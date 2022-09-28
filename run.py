import logging
import os
import sys
import traceback
from multiprocessing import freeze_support
from pathlib import Path
from typing import List, Tuple, Union

from src.config_handling import command_line_tools as cli, parse_args
from src.config_handling.formatting import PrimerParams
from src.config_handling.get_parameters_script import get_config_file, \
    get_config_names
from src.config_handling.length_compatability import potential_length_overflow
from src.hetero_spacer_generator.primer_tools import hms
from src.multiplex_spacer_generator.binding_align import get_time_string
from src.multiplex_spacer_generator.find_best_spacer_combo import \
    FindSpacerCombo
from src.multiplex_spacer_generator.generate_spacer_seqs import \
    get_best_heterogeneity_spacer_seqs_threadable
from src.multiplex_spacer_generator.primer_pool import PrimerPool
from src.seq_alignment_analyser.find_binding_pairs import FindBindingPairs
from src.seq_alignment_analyser.iterator_manager import RESTRICTED
from time import time
from datetime import datetime

from src.seq_alignment_analyser.sequence_management import BindingPair, rev_comp

DIV = '\n' + '=' * 80 + '\n'

# Logger config.
log = logging.getLogger('root')
prog_log = logging.Logger('prog_log', level=0)
prog_log.addHandler(logging.StreamHandler(sys.stdout))


def get_secs_left(end_time: Union[float, int]) -> int:
    """Returns the number of seconds until <end_time>"""
    return int(end_time - time())


class RunController:
    """
    Responsible for storing and executing high-level program functionality.
    """

    num_reps: int
    rep_number: int
    warn: bool

    def __init__(self) -> None:
        """Fetch a config file from the user and initialise attributes required
        to execute the run specified in that config."""
        ns = parse_args.get_cla_namespace()
        try:
            config = ns.config[:-3]
        except TypeError:
            config = get_config_file()

        if config not in get_config_names():
            raise ValueError(config + ' is not a valid config file.')

        # Get config file.
        config_file_name = 'configs.' + config
        self.config = __import__(config_file_name, fromlist=[''])
        parse_args.get_modified_config(self.config)
        self.num_reps = self.config.num_repetitions
        self.rep_number = 0

        # Test that outfile parent directory exists.
        if self.config.out_filepath and \
                not Path(self.config.out_filepath).parent.exists():
            raise ValueError(str(Path(self.config.out_filepath).parent) +
                             ' does not exist.')

        # Logger Configuration
        if self.config.silent:
            prog_log.setLevel(40)

        dir_path = Path(os.path.dirname(__file__))
        now = datetime.now()
        dt_string = now.strftime("%b-%d-%Y %H:%M:%S")
        log_filename = config_file_name + ' ' + dt_string
        log.addHandler(logging.FileHandler(dir_path / 'logs' / log_filename))
        if self.config.verbose:
            log.addHandler(logging.StreamHandler(sys.stdout))
        log.info("Starting...")

        # Runtime Management.
        self.start_time = time()
        self.end_time = int(time() + hms(*self.config.runtime_estimate))

        self.write_out_title()

    def reset_time_params(self) -> None:
        """Sets the runtime restrictions."""
        # Runtime Management.
        self.start_time = time()
        self.end_time = int(time() + hms(*self.config.runtime_estimate))

    def get_config_type(self) -> str:
        """Returns the config type of the user-selected config file."""
        return self.config.config_type

    def write_out_title(self) -> None:
        """Writes a title to the outfile."""
        title = '\n' + '=' * 30 + datetime.now().strftime("%b-%d-%Y %H:%M:%S") \
                + '=' * 30 + '\n\n'
        self.write_to_outfile(title)

    def write_to_outfile(self, s: str) -> None:
        """Writes <s> to the outfile if it exists."""
        if self.config.out_filepath:
            if Path(self.config.out_filepath).exists():
                mode = 'a'
            else:
                mode = 'w'

            with open(self.config.out_filepath, mode) as outfile:
                outfile.write(s)

    def get_heterogeneity_spacers(self) -> None:

        f_binding_seqs = [pair.forward for pair in self.config.binding_sequences]
        r_binding_seqs = [pair.reverse for pair in self.config.binding_sequences]
        targ_names = ['Pair #' + str(i + 1)
                      for i in range(len(self.config.binding_sequences))]

        f_combo, r_combo = self._get_spacers(f_binding_seqs, r_binding_seqs)

        f_5p = [adapter_pair[0] for adapter_pair in self.config.adapters]
        r_5p = [adapter_pair[1] for adapter_pair in self.config.adapters]

        binding_set = self.get_heterogeneity_spacer_sequences(
            f_5p, f_binding_seqs, r_5p, r_binding_seqs, f_combo, r_combo,
            get_secs_left(self.end_time))

        binding_set.add_targ_names(targ_names)

        self.write_to_outfile(DIV + self.get_rep_str() + str(binding_set))
        prog_log.info(binding_set)

    def get_rep_str(self) -> str:
        return "Run #" + str(self.rep_number) + '\n'

    def get_binding_regions(self) -> None:
        """

        Finds valid binding regions and returns them to the user.

        """

        self.potentially_too_long_primers_check()

        best_binding_params = self._find_binding_regions()

        s = 'Binding Parameters:\n\n'

        s += binding_selection_str(best_binding_params)

        prog_log.info(s)

        self.write_to_outfile(DIV + self.get_rep_str() + s)

    def full_run(self) -> None:
        """

        Runs the complete primer binding selection and heterogeneity region
        selection using the information stored in the config file.

        """

        self.potentially_too_long_primers_check()

        best_binding_params = self._find_binding_regions()

        binding_param_str = 'Binding Parameters:\n'

        binding_param_str += binding_selection_str(best_binding_params)

        prog_log.info(binding_param_str)

        f_binding_seqs = [pair.get_f_seq() for pair in best_binding_params]
        r_binding_seqs = [pair.get_r_seq() for pair in best_binding_params]
        targ_names = [pair.target_name for pair in best_binding_params]

        f_combo, r_combo = self._get_spacers(f_binding_seqs, r_binding_seqs)

        f_5p = [adapter_pair[0] for adapter_pair in self.config.adapters]
        r_5p = [adapter_pair[1] for adapter_pair in self.config.adapters]

        binding_set = self.get_heterogeneity_spacer_sequences(
            f_5p, f_binding_seqs, r_5p, r_binding_seqs, f_combo, r_combo,
            get_secs_left(self.end_time))

        binding_set.add_targ_names(targ_names)

        primers_str = 'Full Primer Sequences:\n\n' + str(binding_set) + '\n'

        prog_log.info(primers_str)

        self.write_to_outfile(DIV + self.get_rep_str() + binding_param_str +
                              primers_str)

    def _get_spacers(self, f_binding_seqs: List[str],
                     r_binding_seqs: List[str]) -> Tuple[List[int], List[int]]:
        """This function handles the case when a user wants to iteratively find the
        binding combo. If the user is fine with failing, a specific length can be
        specified."""
        prog_log.info('Finding Heterogeneity Spacers...')

        def get_best(binding_seqs: List[str], alloted_time: int) -> List[int]:
            combo = []
            gb_end_time = time() + alloted_time

            success = True
            max_spacer_size = self.config.max_spacer_length
            while success:

                log.info('\n*** Trying to find spacers with max spacer size: ' +
                         str(max_spacer_size) + '.\n')

                fsc = FindSpacerCombo(runtime=get_secs_left(gb_end_time),
                                      max_threads=self.config.num_threads,
                                      seqs=binding_seqs,
                                      hetero_region_size=self.config.hetero_region_len,
                                      max_spacer_size=max_spacer_size)

                try:
                    combo = fsc.run()
                    log.info('\n*** Successfully found spacers! Decreasing max '
                             'spacer size.\n')
                except RuntimeError:
                    log.info(
                        '\n*** Failed to find spacers. Returning spacers with '
                        'lowest max length found so far.\n')
                    success = False

                max_spacer_size -= 1
            if get_secs_left(gb_end_time) < 0:
                raise RuntimeError(
                    'The program ran out of time to find the best '
                    'spacers. Be sure to allot a longer runtime when '
                    'searching for the min spacers.')

            if not combo:
                raise RuntimeError('Failed to find a spacer combo matching the '
                                   'given specifications.')
            return combo

        if self.config.hetero_region_len == self.config.max_spacer_length:
            f_fsc = FindSpacerCombo(runtime=get_secs_left(self.end_time) // 3,
                                    max_threads=self.config.num_threads,
                                    seqs=f_binding_seqs,
                                    hetero_region_size=self.config.hetero_region_len,
                                    max_spacer_size=self.config.max_spacer_length)

            f_combo = f_fsc.run()

            r_fsc = FindSpacerCombo(runtime=get_secs_left(self.end_time) // 2,
                                    max_threads=self.config.num_threads,
                                    seqs=r_binding_seqs,
                                    hetero_region_size=self.config.hetero_region_len,
                                    max_spacer_size=self.config.max_spacer_length)

            r_combo = r_fsc.run()

        else:
            # For both the reverse and forward spacer, continue decreasing max
            # spacer length until failure.
            f_combo, r_combo = \
                get_best(f_binding_seqs, get_secs_left(self.end_time) // 3), \
                get_best(r_binding_seqs, get_secs_left(self.end_time) // 2)

        complete_time_elapsed_msg(self.start_time)
        return f_combo, r_combo

    def potentially_too_long_primers_check(self) -> None:

        if self.config.silent or self.config.no_warn:
            return

        # Check length requirements.
        potentially_too_long_primer = potential_length_overflow(self.config.adapters,
                                                                self.config.max_spacer_length,
                                                                self.config.binding_region_len)
        if potentially_too_long_primer:
            continue_running = cli.yes_no_prompt(
                'This input may produce primers of length ' +
                str(potentially_too_long_primer) +
                ' which exceeds the 60 base pair limit as imposed by '
                'primer3. In order to continue, primer structures will be  '
                'truncated at their 5\' ends down to meet the  length criteria.'
                ' Enter Y to continue.')
            if not continue_running:
                exit(0)

    def _find_binding_regions(self) -> List[BindingPair]:
        """Gets the binding regions using parameters specified in the config
        file."""
        # Assemble Primer Params

        primer_params = PrimerParams(
            primer_primer_distance=self.config.primer_primer_distance,
            primer_target_distance=self.config.primer_target_distance,
            target_region_len=self.config.target_region_len,
            binding_region_len=self.config.binding_region_len,
            ideal_binding_size=self.config.ideal_binding_size,
            max_binding_target_len=self.config.max_binding_target_len)

        # Find the best regions for primer to bind.

        prog_log.info('Finding Binding Pairs...')
        fbp = FindBindingPairs(target_sites=self.config.targets,
                               adapters=self.config.adapters,
                               primer_params=primer_params,
                               alignments_path=self.config.alignments_path,
                               targ_mt=self.config.target_melting_temp,
                               max_mt_deviance=self.config.max_mt_deviance,
                               aln_type=self.config.alignment_type,
                               do_prog_bars=False,
                               mode=RESTRICTED,
                               num_rand_to_choose_from=self.config.how_random)

        best_binding_params = fbp.get_best_binding_pairs()
        best_binding_params.sort(key=lambda a: a.target_name)

        complete_time_elapsed_msg(self.start_time)

        return best_binding_params

    def get_heterogeneity_spacer_sequences(self, f_5p: List[str],
                                           f_binding_seqs: List[str],
                                           r_5p: List[str],
                                           r_binding_seqs: List[str],
                                           f_spacer_combo: List[int],
                                           r_spacer_combo: List[int],
                                           runtime: int)\
            -> PrimerPool:
        """Spends <runtime> looking for a PrimerPool constructed with given
        components and randomly generated valid heterogeneity spacers with
        minimal dGs among some number of the worst structures.."""

        # Average the worst half of structures to calculate the score for a set.
        num_structs_to_avg = len(f_5p) * len(r_5p) // 2

        prog_log.info('Optimizing Heterogeneity Spacer Sequences...')
        binding_set = get_best_heterogeneity_spacer_seqs_threadable(
            f_5p=f_5p,
            f_binding=f_binding_seqs,
            r_5p=r_5p,
            r_binding=r_binding_seqs,
            f_spacers=f_spacer_combo,
            r_spacers=r_spacer_combo,
            allowed_seconds=runtime,
            num_structs_to_avg=num_structs_to_avg,
            num_threads=self.config.num_threads)

        complete_time_elapsed_msg(self.start_time)

        return binding_set


def main():

    # Configure Run Controller
    run_control = RunController()
    config_type = run_control.get_config_type()

    # If doing multiple replicates, store progress.
    num_failures = 0
    num_successes = 0
    err_str = ''

    # Execute specified number of replicates.
    prog_log.info_rep_number = run_control.num_reps > 1
    for r in range(run_control.num_reps):
        run_control.reset_time_params()
        # Log repetition info.
        if prog_log.info_rep_number:
            msg = 'Beginning repetition #' + str(r + 1) + '.'
            log.info(msg)
            prog_log.info(msg)

        run_control.rep_number += 1
        # Don't ask for user for input in later replicates.
        if r > 0:
            run_control.config.no_warn = True

        # Attempt a run. If an error is raised, return it to the user but
        # continue the run.
        try:
            # Execute appropriate run type.
            match config_type:

                case 'full':
                    run_control.full_run()

                case 'hetero':
                    run_control.get_heterogeneity_spacers()

                case 'binding':
                    run_control.get_binding_regions()
            num_successes += 1

        except Exception as e:
            if run_control.num_reps > 1:
                tb = traceback.format_exc()
                prog_log.info(tb)

                num_failures += 1
                err_str += str(e) + '\n'
            else:
                num_failures += 1
                err_str += str(e) + '\n'
                run_control.write_to_outfile(
                    'Num Failures = ' + str(num_failures) +
                    '\n' + 'Failures: \n' + err_str)
                raise e
        if num_successes == run_control.config.max_successes:
            break

    run_control.write_to_outfile('Num Failures = ' + str(num_failures) +
                                 '\n' + 'Failures: \n' + err_str)


def binding_selection_str(binding_params: List[BindingPair]) -> str:
    """Visualises the given binding params. """
    bp_strs = ''

    for bp in binding_params:
        bp_str = ''.join(
            [
                '>', bp.target_name, '\n',
                '\t Conservation Score: ', str(bp.get_conservation_score()),
                '\n',
                'Forward Binding Region: \n',
                '\tBounds: ', str(bp.f_5p), ' - ', str(bp.f_5p + bp.f_len - 1),
                '\n',
                '\tSequence: ', bp.get_f_seq(), '\n',
                'Reverse Binding Region: \n',
                '\tBounds: ', str(bp.r_5p - bp.r_len + 1), ' - ', str(bp.r_5p),
                '\n',
                '\tSequence: ', bp.get_r_seq(), '\n',
                '\tReverse Compliment: ', rev_comp(bp.get_r_seq()), '\n\n'
            ]
        )
        bp_strs += bp_str

    return bp_strs


def complete_time_elapsed_msg(start_time: float) -> None:
    """Prints a message alerting the user to the time elapsed since
    <start_time>."""
    ta = int(time() - start_time)
    prog_log.info('Complete. Time elapsed: \n' + get_time_string(ta))


if __name__ == '__main__':
    freeze_support()
    main()
