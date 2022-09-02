from src.multiplex_spacer_generator.find_best_spacer_combo import FindSpacerCombo
from multiplex_spacer_generator.generate_spacer_seqs import *
import random
from test_files.fixtures_and_helpers import add_stdout_handler, \
    configure_log_out
from test_files.multiplex_testing.test_generate_spacer_seqs import rand_seq

NUM_THREADS = 8

add_stdout_handler()

configure_log_out('stress_test.py')


def test_hetero_seq_generation_threaded() -> None:
    max_het = 12
    max_seqs = 12
    num_reps = 20
    max_adapter_len = 24
    max_seq_len = 24
    spacer_finder_runtime = 10
    seq_generation_runtime = 10

    for _ in range(num_reps):
        het_len = random.randint(6, max_het)
        num_seqs = random.randint(6, max_seqs)
        f_binding = [rand_seq(max_seq_len) for _ in range(num_seqs)]
        r_binding = [rand_seq(max_seq_len) for _ in range(num_seqs)]
        f_adapters = [rand_seq(max_adapter_len) for _ in range(num_seqs)]
        r_adapters = [rand_seq(max_adapter_len) for _ in range(num_seqs)]
        num_structs_to_view = num_seqs ** 2 // 2

        spacer_finder = FindSpacerCombo(spacer_finder_runtime, NUM_THREADS,
                                        f_binding, het_len)
        f_spacers = spacer_finder.run()
        spacer_finder = FindSpacerCombo(spacer_finder_runtime, NUM_THREADS,
                                        r_binding, het_len)
        r_spacers = spacer_finder.run()

        final_pool = get_best_heterogeneity_spacer_seqs_threadable(
            f_adapters, f_binding, r_adapters, r_binding, f_spacers, r_spacers,
            seq_generation_runtime, num_structs_to_view, NUM_THREADS)
        assert final_pool
