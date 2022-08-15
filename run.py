import logging
import os
import sys
from collections.abc import Sequence
from multiprocessing import freeze_support
from pathlib import Path

from config_handling.formatting import PrimerParams
from config_handling.get_parameters_script import get_config_file
from hetero_spacer_generator.primer_tools import hms
from multiplex_spacer_generator.find_best_spacer_combo import FindSpacerCombo, \
    MIN_PER_PROCESS
from multiplex_spacer_generator.generate_spacer_seqs import \
    get_best_heterogeneity_spacer_seqs_threadable
from multiplex_spacer_generator.primer_pool import PrimerPool
from seq_alignment_analyser.find_binding_pairs import FindBindingPairs
from seq_alignment_analyser.iterator_manager import RESTRICTED
from time import time
from datetime import datetime

def main():
    config_file_name = 'configs.' + get_config_file()
    config = __import__(config_file_name, fromlist=[''])

    # Logger config.

    log = logging.getLogger('root')
    log.setLevel(0)
    dir_path = Path(os.path.dirname(__file__))
    now = datetime.now()
    dt_string = now.strftime("%b-%d-%Y %H:%M:%S")
    log_filename = config_file_name[:-3] + ' ' + dt_string
    log.addHandler(logging.FileHandler(dir_path / log_filename))
    log.addHandler(logging.StreamHandler(sys.stdout))
    log.info("Starting...")

    # Runtime Management.

    end_time = int(time() + hms(*config.runtime_estimate))
    get_secs_left = lambda : int(end_time - time())

    # Assemble Primer Params

    primer_params = PrimerParams(primer_primer_distance=config.primer_primer_distance,
                                 primer_target_distance=config.primer_target_distance,
                                 target_region_len=config.target_region_len,
                                 binding_region_len=config.binding_region_len,
                                 ideal_binding_size=config.ideal_binding_size,
                                 max_binding_target_len=config.max_binding_target_len)


    # Find the best regions for primer to bind.

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
    f_binding_seqs = [pair.get_f_seq() for pair in best_binding_params]
    r_binding_seqs = [pair.get_r_seq() for pair in best_binding_params]
    targ_names = [pair.target_name for pair in best_binding_params]

    # Find the best spacer lengths that produce the requested level of
    # heterogeneity.

    f_fsc = FindSpacerCombo(runtime=get_secs_left() // 3,
                          max_threads=config.num_threads,
                          seqs=f_binding_seqs,
                          hetero_region_size=config.hetero_region_len)

    f_combo = f_fsc.run()

    r_fsc = FindSpacerCombo(runtime=get_secs_left() // 2,
                          max_threads=config.num_threads,
                          seqs=r_binding_seqs,
                          hetero_region_size=config.hetero_region_len)

    r_combo = r_fsc.run()


    # Find the best sequences to insert into the heterogeneity regions in order to
    # reduce dimerisation potential.


    f_5p = [adapter_pair[0] for adapter_pair in config.adapters]
    r_5p = [adapter_pair[1] for adapter_pair in config.adapters]

    # Average the worst half of structures to calculate the score for a set.
    num_structs_to_avg = len(f_5p) * len(r_5p) // 2

    binding_set = get_best_heterogeneity_spacer_seqs_threadable(
                                    f_5p=f_5p,
                                    f_binding=f_binding_seqs,
                                    r_5p=r_5p,
                                    r_binding=r_binding_seqs,
                                    f_spacers=f_combo,
                                    r_spacers=r_combo,
                                    allowed_seconds=get_secs_left(),
                                    num_structs_to_avg=num_structs_to_avg,
                                    num_threads=config.num_threads)

    binding_set.add_targ_names(targ_names)
    print(binding_set)

if __name__ == '__main__':
    freeze_support()
    main()