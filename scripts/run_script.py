import os
import sys

config = sys.argv[1]
out_file = sys.argv[2]

how_random = 10
num_repetitions = 1
max_successes = 10
runtime_estimate = '0 5 0'
num_threads = 8

out_filepath = ''
hetero_region_length = 0
max_spacer_length = 0


def run_cmd():
    cmd = [
        'python3.10 ', 'run.py',
        '--config', config,
        '--how_random', how_random,
        '--num_repetitions', num_repetitions,
        '--runtime_estimate', runtime_estimate,
        # '--verbose',
        '--no_warn',
        # '--silent',
        '--hetero_region_length', hetero_region_length,
        '--num_threads', num_threads,
        '--max_spacer_length', max_spacer_length,
        '--out_filepath', out_filepath,
        '--max_successes', max_successes
    ]

    cmd = [str(s) for s in cmd]
    cmd = ' '.join(cmd)

    print('Executing: ', cmd)

    os.system(cmd)


for total_h_len in range(8, 12):
    for max_het_spacer in range(6, 10):
        hetero_region_length = total_h_len
        max_spacer_length = max_het_spacer

        out_filepath = out_file + '/' + str(max_spacer_length) + '_' + \
                       str(hetero_region_length) + '.txt'

        run_cmd()
