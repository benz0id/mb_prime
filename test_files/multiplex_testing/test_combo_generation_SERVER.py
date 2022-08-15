import os

from hetero_spacer_generator.primer_tools import hms
from test_files.fixtures_and_helpers import configure_log_out, gen_random_seq, \
    TEST_PATH, \
    TEST_OUTPUT_PATH, add_stdout_handler
from multiplex_spacer_generator.find_best_spacer_combo import *
import pathlib
import pytest

add_stdout_handler()

MAX_NUM_PROCS = 8
TOTAL_RUNTIME = hms(15, 0, 0)

exp_time = time() + TOTAL_RUNTIME

rem_time = lambda: int(exp_time - time())

log_num = 0
def new_log():
    global log_num
    configure_log_out('FindSpacerCombo_' + str(log_num))
    log_num += 1

def test_run_full_few_seqs() -> None:
    configure_log_out('test_run_full_many_seqs')
    basic_fsc = FindSpacerCombo(rem_time(), MAX_NUM_PROCS,
                                [
                                    'ATGCATGCATGC',
                                    'CATGCATGCATGC',
                                    'GCATGCATGCATGC',
                                    'TGCATGCATGCATGC',
                                ],
                                12)
    combo = basic_fsc.run()
    assert combo


def test_run_full_many_random_seqs_len_9() -> None:
    num_seqs = 9
    hetero = 12
    num_iter = 10
    for i in range(num_iter):
        configure_log_out(str(hetero) + '^' + str(num_seqs) + '#' + str(i + 1))
        basic_fsc = FindSpacerCombo(hms(0, 5, 0), MAX_NUM_PROCS,
                                    [str(gen_random_seq(hetero))
                                     for _ in range(num_seqs)],
                                    12)
        combo = basic_fsc.run()
        assert combo


def test_run_full_many_random_seqs_len_8() -> None:
    num_seqs = 8
    hetero = 12
    num_iter = 10
    for i in range(num_iter):
        configure_log_out(str(hetero) + '^' + str(num_seqs) + '#' + str(i + 1))
        basic_fsc = FindSpacerCombo(hms(0, 5, 0), MAX_NUM_PROCS,
                                    [str(gen_random_seq(hetero))
                                     for _ in range(num_seqs)],
                                    12)
        combo = basic_fsc.run()
        assert combo


def stress_test() -> None:
    num_seqs = 20
    hetero = 20
    num_iter = 10
    for i in range(num_iter):
        configure_log_out(str(hetero) + '^' + str(num_seqs) + '#' + str(i + 1))
        basic_fsc = FindSpacerCombo(rem_time(), MAX_NUM_PROCS,
                                    [str(gen_random_seq(hetero))
                                     for _ in range(num_seqs)],
                                    12)
        combo = basic_fsc.run()
        assert combo


if __name__ == '__main__':
    cur_path = __file__
    print(cur_path)
    cmd = ''.join(
        [
            'pytest ', str(cur_path),
            ' --log-cli-level=10  '
            # ' > ', str(TEST_OUTPUT_PATH / 'pytest_summary.log')
        ])
    print(cmd)
    os.system(cmd)
