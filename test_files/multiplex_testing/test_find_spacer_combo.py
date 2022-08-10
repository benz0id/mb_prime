from test_files.fixtures_and_helpers import configure_log_out, gen_random_seq
from multiplex_spacer_generator.find_best_spacer_combo import *

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

    seqs = [
        'ATGCATGCATGC',
        'CATGCATGCATGC',
        'GCATGCATGCATGC',
        'TGCATGCATGCATGC',
        'ATGCATGCATGC',
        'CATGCATGCATGC',
        'GCATGCATGCATGC',
        'TGCATGCATGCATGC'
    ]
    basic_fsc = FindSpacerCombo(1800, 8, seqs, 12)
    combo = basic_fsc.run()
    combo_str(combo, seqs)
    assert combo
