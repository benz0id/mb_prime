from test_files.fixtures_and_helpers import configure_log_out, gen_random_seq
from multiplex_spacer_generator.find_best_spacer_combo import *
import random

from test_files.multiplex_testing.test_find_binding_pairs import rand_seq

configure_log_out('FindSpacerCombo')


def test_thread_params() -> None:
    basic_fsc = FindSpacerCombo(60, 5, ['ATCG', 'GATC', 'CCCC'], 3, min_per_process=2)
    thread_params = basic_fsc._get_thread_params(3)
    assert thread_params == [
        (2, [3, 0, 0]),
        (2, [1, 2, 0]),
        (2, [2, 0, 1]),
        (2, [0, 2, 1]),
        (2, [0, 1, 2])
    ]


def test_thread_params_upper() -> None:
    basic_fsc = FindSpacerCombo(60, 5, ['ATCG', 'GATC', 'CCCC'], 4,
                                min_per_process=2)
    thread_params = basic_fsc._get_thread_params(12)
    assert thread_params == \
           [
               (1, [4, 4, 4])
           ]


def test_thread_params_lower() -> None:
    basic_fsc = FindSpacerCombo(60, 5, ['ATCG', 'GATC', 'CCCC'], 4,
                                min_per_process=2)
    thread_params = basic_fsc._get_thread_params(0)
    assert thread_params == \
           [
               (1, [0, 0, 0])
           ]


def test_thread_param_complex() -> None:
    basic_fsc = FindSpacerCombo(1000000, 10,
                                [str(gen_random_seq(12)) for _ in range(8)],
                                3, min_per_process=1)
    thread_params = basic_fsc._get_thread_params(1)
    assert len(thread_params) == 8

    basic_fsc = FindSpacerCombo(1000000, 28,
                                [str(gen_random_seq(12)) for _ in range(8)],
                                3, min_per_process=1)
    thread_params = basic_fsc._get_thread_params(2)

    assert len(thread_params) == 28


def test_run_full_single_proc() -> None:
    basic_fsc = FindSpacerCombo(20, 2,
                                [
                                    'ATGCATGCATGC',
                                    'CATGCATGCATGC',
                                    'GCATGCATGCATGC',
                                    'TGCATGCATGCATGC'
                                ],
                                12, min_per_process=1)
    combo = basic_fsc.run()
    assert combo == [0, 0, 0, 0]

def test_run_full_many_seqs_single_proc() -> None:
    basic_fsc = FindSpacerCombo(10, 1,
                                [
                                    'ATGCATGCATGC',
                                    'CATGCATGCATGC',
                                    'GCATGCATGCATGC',
                                    'TGCATGCATGCATGC',
                                    'ATGCATGCATGC',
                                    'CATGCATGCATGC',
                                    'GCATGCATGCATGC',
                                    'TGCATGCATGCATGC'
                                ],
                                12)
    combo = basic_fsc.run()
    assert combo


def test_run_full_many_seqs_multi_proc() -> None:

    seqs = [str(gen_random_seq(12)) for _ in range(8)]
    basic_fsc = FindSpacerCombo(1800, 8, seqs, 12)
    combo = basic_fsc.run()
    print('\n', combo_str(combo, seqs))
    assert combo


def test_binding_align_iteration() -> None:
    max_het = 20
    max_seqs = 5
    num_reps = 100
    max_seq_len = 24
    num_incr = 1000

    for _ in range(num_reps):
        het_len = random.randint(0, max_het)
        num_seqs = random.randint(1, max_seqs)
        f_binding = [rand_seq(max_seq_len) for _ in range(num_seqs)]
        spacer_sizes = [random.randint(0, het_len) for _ in range(num_seqs)]
        binding_align = NumpyBindingAlign(f_binding, spacer_sizes, het_len)

        for _ in range(num_incr):
            binding_align.incr_spacer_sizes()
            assert min(binding_align._spacer_sizes) >= 0
            assert max(binding_align._spacer_sizes) <= binding_align._num_hetero


