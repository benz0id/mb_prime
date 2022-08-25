import logging
import os
import sys
from collections.abc import Sequence
from multiprocessing import freeze_support
from pathlib import Path
from typing import List, Tuple

import config_handling.command_line_tools as cli
from config_handling.formatting import PrimerParams
from config_handling.get_parameters_script import get_config_file
from config_handling.length_compatability import potential_length_overflow
from hetero_spacer_generator.primer_tools import hms
from multiplex_spacer_generator.binding_align import SpacerCombo, \
    get_time_string
from multiplex_spacer_generator.find_best_spacer_combo import FindSpacerCombo, \
    MIN_PER_PROCESS
from multiplex_spacer_generator.generate_spacer_seqs import \
    get_best_heterogeneity_spacer_seqs_threadable
from multiplex_spacer_generator.primer_pool import PrimerPool
from seq_alignment_analyser.find_binding_pairs import FindBindingPairs
from seq_alignment_analyser.iterator_manager import RESTRICTED
from time import time
from datetime import datetime

from seq_alignment_analyser.sequence_management import BindingPair, rev_comp

get_secs_left = lambda end_time: int(end_time - time())


def main():
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

    # Check length requirements.
    potentially_too_long_primer = potential_length_overflow(config.adapters,
                                                            config.max_spacer_length,
                                                            config.binding_region_len)
    if potentially_too_long_primer:
        continue_running = cli.yes_no_prompt(
            'This input may produce primers of length ' +
            str(potentially_too_long_primer) + ' which exceeds the 60 base pair'
                                               ' limit as imposed by primer3. In order to continue, primer '
                                               'structures will be  truncated at their 5\' ends down to meet the '
                                               'length criteria. Enter Y to continue.')
        if not continue_running:
            exit(0)

    # Runtime Management.

    start_time = time()
    end_time = int(time() + hms(*config.runtime_estimate))

    # Assemble Primer Params

    primer_params = PrimerParams(primer_primer_distance=config.primer_primer_distance,
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

    f_binding_seqs = [pair.get_f_seq() for pair in best_binding_params]
    r_binding_seqs = [pair.get_r_seq() for pair in best_binding_params]
    targ_names = [pair.target_name for pair in best_binding_params]

    # Find the best spacer lengths that produce the requested level of
    # heterogeneity.

    prog_log.info('Finding Heterogeneity Spacers...')
    f_combo, r_combo = get_spacers(config_file_name, end_time, f_binding_seqs,
                                   r_binding_seqs, log)
    complete_time_elapsed_msg(start_time)


    # Find the best sequences to insert into the heterogeneity regions in order to
    # reduce dimerisation potential.


    f_5p = [adapter_pair[0] for adapter_pair in config.adapters]
    r_5p = [adapter_pair[1] for adapter_pair in config.adapters]

    # Average the worst half of structures to calculate the score for a set.
    num_structs_to_avg = len(f_5p) * len(r_5p) // 2

    prog_log.info('Optimizing Heterogeneity Spacer Sequences...')
    binding_set = get_best_heterogeneity_spacer_seqs_threadable(
                                    f_5p=f_5p,
                                    f_binding=f_binding_seqs,
                                    r_5p=r_5p,
                                    r_binding=r_binding_seqs,
                                    f_spacers=f_combo,
                                    r_spacers=r_combo,
                                    allowed_seconds=get_secs_left(end_time),
                                    num_structs_to_avg=num_structs_to_avg,
                                    num_threads=config.num_threads)

    binding_set.add_targ_names(targ_names)

    print(binding_selection_str(best_binding_params))

    print(binding_set)


def binding_selection_str(binding_params: List[BindingPair]) -> str:
    """Visualises the given binding params. """
    bp_strs = ''

    for bp in binding_params:
        bp_str = ''.join(
            [
                '>', bp.target_name, '\n'
                'Forward Binding Region: \n',
                '\t Bounds: ', str(bp.f_5p), ' - ', str(bp.f_5p + bp.f_len - 1),
                '\n',
                '\t Sequence: ', bp.get_f_seq(), '\n',
                'Reverse Binding Region: \n',
                '\t Bounds: ', str(bp.r_5p), ' - ', str(bp.r_5p + bp.r_len - 1),
                '\n',
                '\t Sequence: ', bp.get_r_seq(), '\n',
                '\t Reverse Compliment: ', rev_comp(bp.get_r_seq()), '\n\n'
            ]
        )
        bp_strs += bp_str

    return bp_strs


def complete_time_elapsed_msg(start_time: float) -> None:
    """Prints a message alerting the user to the time elapsed since
    <start_time>."""
    ta = int(time() - start_time)
    print('Complete. Time elapsed: ', get_time_string(ta))


def get_spacers(config_file_name: str, end_time: int,
                f_binding_seqs: List[str], r_binding_seqs: List[str],
                log: logging.Logger) \
    -> Tuple[List[int], List[int]]:
    """This function handles the case when a user wants to iteratively find the
    binding combo. If the user is fine with failing, a specific length can be
    specified."""
    config = __import__(config_file_name, fromlist=[''])

    def get_best(binding_seqs: List[str], alloted_time: int) -> List[int]:
        combo = []
        gb_end_time = time() + alloted_time

        success = True
        max_spacer_size = config.hetero_region_len
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
                log.info('\n*** Failed to find spacers. Returning spacers with '
                         'lowest max length found so far.\n')
                success = False

            max_spacer_size -= 1
        if get_secs_left(gb_end_time) < 0:
            raise RuntimeError('The program ran out of time to find the best '
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
    return f_combo, r_combo


if __name__ == '__main__':
    freeze_support()
    main()
