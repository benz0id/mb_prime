import logging
import os
import sys
from multiprocessing import freeze_support
from pathlib import Path
from typing import List, Tuple, Union

from src.config_handling import command_line_tools as cli
from src.config_handling.formatting import PrimerParams
from src.config_handling.get_parameters_script import get_config_file
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

# Freeze support check. Avoids recursive program execution when multithreading.
if __name__ == '__main__':
    freeze_support()

# Get config file.
config_file_name = 'configs.' + get_config_file()
config = __import__(config_file_name, fromlist=[''])

# Logger config.
prog_log = logging.getLogger('prog_log')
prog_log.addHandler(logging.StreamHandler(sys.stdout))
prog_log.setLevel(0)

log = logging.getLogger('root')
log.setLevel(0)
dir_path = Path(os.path.dirname(__file__))
now = datetime.now()
dt_string = now.strftime("%b-%d-%Y %H:%M:%S")
log_filename = config_file_name + ' ' + dt_string
log.addHandler(logging.FileHandler(dir_path / 'logs' / log_filename))
if config.verbose:
    log.addHandler(logging.StreamHandler(sys.stdout))
log.info("Starting...")

# Runtime Management.
start_time = time()
end_time = int(time() + hms(*config.runtime_estimate))


def get_secs_left(end_time: Union[float, int]) -> int:
    """Returns the number of seconds until <end_time>"""
    return int(end_time - time())


class RunController:
    """





    """

    def get_heterogeneity_spacers(self) -> None:
        f_binding_seqs = [pair.get_f_seq() for pair in ]
        r_binding_seqs = [pair.get_r_seq() for pair in ]
        targ_names = [pair.target_name for pair in ]

        f_combo, r_combo = self._get_spacers(f_binding_seqs, r_binding_seqs)

        f_5p = [adapter_pair[0] for adapter_pair in config.adapters]
        r_5p = [adapter_pair[1] for adapter_pair in config.adapters]

        binding_set = self.get_heterogeneity_spacer_sequences(
            f_5p, f_binding_seqs, r_5p, r_binding_seqs, f_combo, r_combo,
            get_secs_left(end_time))

        binding_set.add_targ_names(targ_names)

    def get_binding_regions(self) -> None:
        """

        Finds valid binding regions and returns them to the user.

        """

        potentially_too_long_primers_check()

        best_binding_params = self._find_binding_regions()

        print('Binding Parameters:\n')

        print(binding_selection_str(best_binding_params))

    def full_run(self) -> None:
        """

        Runs the complete primer binding selection and heterogeneity region
        selection using the information stored in the config file.

        """

        potentially_too_long_primers_check()

        best_binding_params = self._find_binding_regions()

        f_binding_seqs = [pair.get_f_seq() for pair in best_binding_params]
        r_binding_seqs = [pair.get_r_seq() for pair in best_binding_params]
        targ_names = [pair.target_name for pair in best_binding_params]

        f_combo, r_combo = self._get_spacers(f_binding_seqs, r_binding_seqs)

        f_5p = [adapter_pair[0] for adapter_pair in config.adapters]
        r_5p = [adapter_pair[1] for adapter_pair in config.adapters]

        binding_set = self.get_heterogeneity_spacer_sequences(
            f_5p, f_binding_seqs, r_5p, r_binding_seqs, f_combo, r_combo,
            get_secs_left(end_time))

        binding_set.add_targ_names(targ_names)

        print('Binding Parameters:\n')

        print(binding_selection_str(best_binding_params))

        print('Full Primer Sequences:\n')

        print(binding_set)

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
            max_spacer_size = config.max_spacer_length
            while success:

                log.info('\n*** Trying to find spacers with max spacer size: ' +
                         str(max_spacer_size) + '.\n')

                fsc = FindSpacerCombo(runtime=get_secs_left(gb_end_time),
                                      max_threads=config.num_threads,
                                      seqs=binding_seqs,
                                      hetero_region_size=config.hetero_region_len,
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
                    'spacers. Be sure to alot a longer runtime when '
                    'searching for the min spacers.')

            assert combo
            return combo

        if not config.min_spacer_length:
            f_fsc = FindSpacerCombo(runtime=get_secs_left(end_time) // 3,
                                    max_threads=config.num_threads,
                                    seqs=f_binding_seqs,
                                    hetero_region_size=config.hetero_region_len,
                                    max_spacer_size=config.max_spacer_length)

            f_combo = f_fsc.run()

            r_fsc = FindSpacerCombo(runtime=get_secs_left(end_time) // 2,
                                    max_threads=config.num_threads,
                                    seqs=r_binding_seqs,
                                    hetero_region_size=config.hetero_region_len,
                                    max_spacer_size=config.max_spacer_length)

            r_combo = r_fsc.run()

        else:
            # For both the reverse and forward spacer, continue decreasing max
            # spacer length until failure.
            f_combo, r_combo = \
                get_best(f_binding_seqs, get_secs_left(end_time) // 3), \
                get_best(r_binding_seqs, get_secs_left(end_time) // 2)

        complete_time_elapsed_msg(start_time)
        return f_combo, r_combo

    def _find_binding_regions(self) -> List[BindingPair]:
        """Gets the binding regions using parameters specified in the config
        file."""
        # Assemble Primer Params

        primer_params = PrimerParams(
            primer_primer_distance=config.primer_primer_distance,
            primer_target_distance=config.primer_target_distance,
            target_region_len=config.target_region_len,
            binding_region_len=config.binding_region_len,
            ideal_binding_size=config.ideal_binding_size,
            max_binding_target_len=config.max_binding_target_len)

        # Find the best regions for primer to bind.

        prog_log.info('Finding Binding Pairs...')
        fbp = FindBindingPairs(target_sites=config.target_sites,
                               adapters=config.adapters,
                               primer_params=primer_params,
                               alignments_path=config.alignments_path,
                               targ_mt=config.target_melting_temp,
                               max_mt_deviance=config.max_mt_deviance,
                               aln_type=config.alignment_type,
                               do_prog_bars=False,
                               mode=RESTRICTED)

        best_binding_params = fbp.get_best_binding_pairs()
        best_binding_params.sort(key=lambda a: a.target_name)

        complete_time_elapsed_msg(start_time)

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
            num_threads=config.num_threads)

        complete_time_elapsed_msg(start_time)

        return binding_set


def main():
    run_contol = RunController()

    match config.config_type:

        case 'full':
            run_contol.full_run()

        case 'hetero':
            run_contol.get_heterogeneity_spacers()

        case 'binding':
            run_contol._find_binding_regions()






def potentially_too_long_primers_check() -> None:
    # Check length requirements.
    potentially_too_long_primer = potential_length_overflow(config.adapters,
                                                            config.max_spacer_length,
                                                            config.binding_region_len)
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
    print('Complete. Time elapsed: ', get_time_string(ta))


if __name__ == '__main__':
    main()
