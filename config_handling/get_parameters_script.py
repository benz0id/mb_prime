import re
from typing import List
import config_handling.command_line_tools as cli
import config_handling.input_validator as iv
import sys
import config_handling.parameters as param
import os
from pathlib import Path
from seq_alignment_analyser.align import MSA
from test_files.fixtures_and_helpers import ALIGNMENTS_PATH

CONFIG_PATH = Path(os.path.dirname(__file__)).parent / 'configs'

HEADER = ('from config_handling.formatting import *\n'
          'from pathlib import Path\n'
          '\n'
          'verbose = False\n'
          '\n'
          'min_spacer_length = False\n')

AUTOFILL = True


def get_config_file() -> str:
    """Lets the user select an existing config, or creates a new one."""
    configs = [str(filename) for filename in os.listdir(CONFIG_PATH)]
    if '__pycache__' in configs:
        configs.remove('__pycache__')
    configs = [config[:-3] for config in configs]
    configs.insert(0, 'Create New Config')
    sel = cli.menu(
        configs, 'Config File Select', 'Enter the number corresponding to'
                                       'the config file you\'d like to '
                                       'use, or enter 1 to create a new '
                                       'config.')
    if sel == 0:
        return get_new_config()
    else:
        return configs[sel]


def get_new_config() -> str:
    """Prompts the user for the parameters required to run the program. Prints
    these parameters as a series of strings to config file and interpreted as
    python variables. Returns path of selected config file."""
    # Collecting basic configuration variables.

    cli.print_title('Basic Parameters')

    existing_configs = os.listdir(CONFIG_PATH)
    # Remove .py extension
    existing_configs = [name[:-3] for name in existing_configs]

    config_name = param.StrParam(
        'config_name',
        'Enter the name of the config file to be created.',
        [iv.Validation(lambda s: s not in existing_configs, 'Config name '
                                                            'already in use.',
                       is_warning=True),
         iv.Validation(lambda s: not bool(re.findall(r'[^A-Za-z0-9_\-\\]', s)),
                       'Invalid filename')])

    config_out = CONFIG_PATH / (config_name.data + '.py')

    default_ali_path = Path(os.path.dirname(__file__)).parent / 'alignments'
    alignments_path = param.PathParam(
        'alignments_path', 'Enter path to file containing alignments or '
                           '\'DIR\' to use local alignments file.', [],
                            default=default_ali_path)


    alignment_type = param.StrParam('alignment_type',
                                    'Enter the type of alignments'
                                    ' contained in the given'
                                    ' folder',
                                    [iv.VALID_MSA_TYPE])

    target_region_len = param.RangeParam('target_region_len',
    'Enter the allowable lengths of the target region (length of amplicon minus primers).',
    [iv.Validation(lambda s: iv.all_in_range(s, -1, 1000000000),
    'Please enter reasonable lengths (0 <= d).')])

    binding_region_len = param.RangeParam('binding_region_len',
    'Enter the allowable lengths of the binding regions on each primer (number '
    'of bases that bind to target). ',
    [iv.Validation(lambda s: iv.all_in_range(s, -1, 1000000000),
    'Please enter reasonable lengths (0 <= d).')])

    ideal_binding_size = param.IntParam(
        "ideal_binding_size",
    'Enter the ideal length of a binding sequence. Adjustments will be made in'
    ' order to maximise conservation and reach target melting temp.',
    [iv.Validation(lambda s: iv.all_in_range(s, 0, 1000000000),
    'Please enter a reasonable length (0 <= d).')])

    max_binding_target_len = param.IntParam(
        'max_binding_target_len',
        'Enter the maximum combined length of both binding regions and the '
        'target region.',
        [iv.Validation(lambda s: iv.all_in_range(s, 0, 1000000000),
                       'Please enter a reasonable length (0 <= d).')])

    hetero_region_len = param.IntParam(
        'hetero_region_len',
        'Enter the number of bases to be include in the heterogeneity region.',
        [iv.Validation(lambda s: iv.all_in_range(s, 0, binding_region_len.data.stop),
                       'Please enter a reasonable length (0 <= d <= max_binding_region_len).')])

    max_spacer_length = param.IntParam(
        'max_spacer_length',
        'Enter the maximum length of any spacer to be included in a primer in '
        'order to achieve heterogeneity. Setting this too short will often '
        'result in a failure to find valid spacers.',
        [iv.Validation(
            lambda s: iv.all_in_range(s, 0, hetero_region_len.data),
            'Please enter a reasonable length (0 <= d <= hetero_region_len).')])

    # Extract variables required for target region selection from the above
    # params.

    # Twice the length of the smallest allowable primer.
    min_space_between_targets = binding_region_len.data[0] * 2

    min_space_from_ends_of_alignment = binding_region_len.data[0]

    max_target_region_len = target_region_len.data[1]

    # Build MSAs using given parameters.

    name_to_msa = {}
    msa_to_target = {}
    msas = []

    align_filenames = os.listdir(alignments_path.data)

    try:
        for filename in align_filenames:
            align_path = alignments_path.data / filename
            msa = MSA(align_path, filetype=alignment_type.data)

            msas.append(msa)
            msa_to_target[msa] = []
            name_to_msa[filename] = msa
    except IOError:
        cli.eprint(
            'Failed to parse one or more alignments. Ensure that all items'
            'in the given file are valid alignments of the given type. See'
            'logs for more detailed description of failure.')
        quit()

    cli.print_title('Target Selection')

    # Prompt user for targets.
    num_targets = param.IntParam('num_targets',
                                 'Enter the number of targets you'
                                 ' would like to create.', [])

    targets = []

    for target_num in range(1, 1 + num_targets.data):
        sel = cli.menu(align_filenames,
                       'Select Target Alignment',
                       'Enter the number corresponding to the alignment you\'d '
                       'like to target.')
        sel_ali_name = align_filenames[sel]

        new_target = \
            param.TargetRegionParam(sel_ali_name,
                                    max_target_len=max_target_region_len,
                                    min_sep_distance=min_space_between_targets,
                                    other_targets=msa_to_target,
                                    target_number=target_num,
                                    start_end_distance=min_space_from_ends_of_alignment)

        new_target.query_user()

        targets.append(new_target)

    # Prompt user for sequences to be placed 5' of the binding sequences.

    cli.print_title('5\' Primer Sequences')

    print('Enter sequences to be placed 5\' of the binding sequences.')
    adapters = []

    for target_num in range(num_targets.data):
        seq_name = str(targets[target_num].name)
        msg = ''.join(['Enter 5\' sequences for ', seq_name, '.'])
        adapter_pair = param.AdapterParam(seq_name + '_adapters', msg, [])
        adapters.append(adapter_pair)

    cli.print_title('Run Parameter Configuration')

    num_threads = param.IntParam(
        'num_threads',
        'Enter the number of threads you\'d like to use during any '
        'parallelizable stage of the program',
        [iv.Validation(lambda s: 0 < int(s) < os.cpu_count() * 3,
                       'Unreasonable number of Threads requested.')])

    runtime_estimate = param.TimeParam(
        'runtime_estimate', 'Enter a rough estimate of the desired program '
                            'runtime.', [])

    target_melting_temp = param.FloatParam('target_melting_temp',
        'Enter the target melting temp for primers generated by this program.',
        [iv.Validation(lambda s: 20 < float(s) < 100,
        'Please enter a reasonable melting temperature (20 < Tm < 100).')])

    max_mt_deviance = param.FloatParam('max_mt_deviance',
        'Enter the maximum melting temp deviance between primers generated by '
        'this program.',
        [iv.Validation(lambda s: 0 <= float(s) < 40,
        'Please enter a reasonable melting temperature distance (0 <= Tm < 40).')])

    primer_primer_distance = param.IntParam('primer_primer_distance',
        'Enter the minumum distance between any two primer\'s binding sites',
        [iv.Validation(lambda s: 0 <= float(s),
        'Please enter a reasonable distance (0 <= d).')]
                                           )

    primer_target_distance = param.IntParam('primer_target_distance',
    'Enter the minumum distance any primer\'s binding site and any target site.',
    [iv.Validation(lambda s: 0 <= float(s),
    'Please enter a reasonable distance (0 <= d).')])


    # Output the given parameters.

    basic_params = [
        alignments_path, alignment_type, target_region_len,
        binding_region_len, ideal_binding_size, max_binding_target_len,
        hetero_region_len, max_spacer_length, target_melting_temp,
        max_mt_deviance, primer_primer_distance, primer_target_distance,
        runtime_estimate, num_threads
    ]

    out_str = HEADER

    for par in basic_params:
        out_str += (repr(par.get_python_variable_string())[1:-1]) + '\n\n'

    out_str += param.format_as_pylist(targets, 'target_sites') + '\n\n'

    out_str += param.format_as_pylist(adapters, 'adapters') + '\n\n'

    with open(config_out, 'w') as outfile:
        outfile.write(out_str)

    return config_name.data


# Input some predetermined set of strings.
if __name__ == '__main__':
    if not AUTOFILL:
        get_new_config()
        quit()


    class DummyIO:
        i: int
        lines: List[str]

        def __init__(self, lines: List[str]):
            self.i = -1
            self.lines = lines
            self.stdin = sys.stdin
            sys.stdin = self

        def readline(self) -> str:
            self.i += 1
            if self.i == len(self.lines):
                sys.stdin = self.stdin
                return sys.stdin.readline()
            print(self.lines[self.i])
            return self.lines[self.i]


    # The lines to be written to the console automatically.
    to_write = [
        'complex_test_config',
        'Y',
        'DIR',
        'fasta',
        '0, 150',
        '15, 25',
        '20',
        '126',
        '12',
        '8',

        '8',
        '1',
        'targ_1',
        '50, 60',
        '1',
        'targ_2',
        '130, 140',
        '1',
        'targ_3',
        '210, 220',
        '1',
        'targ_4',
        '290, 300',
        '1',
        'targ_5',
        '370, 380',
        '1',
        'targ_6',
        '450, 460',
        '1',
        'targ_7',
        '600, 630',
        '1',
        'targ_8',
        '710, 730',

        'ACACTCTTTCCCTACACGACGCTCTTCCGATCT',
        'GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT',
        'ACACTCTTTCCCTACACGACGCTCTTCCGATCT',
        'GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT',
        'ACACTCTTTCCCTACACGACGCTCTTCCGATCT',
        'GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT',
        'ACACTCTTTCCCTACACGACGCTCTTCCGATCT',
        'GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT',
        'ACACTCTTTCCCTACACGACGCTCTTCCGATCT',
        'GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT',
        'ACACTCTTTCCCTACACGACGCTCTTCCGATCT',
        'GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT',
        'ACACTCTTTCCCTACACGACGCTCTTCCGATCT',
        'GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT',
        'ACACTCTTTCCCTACACGACGCTCTTCCGATCT',
        'GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT',

        '8',
        '0:5:0',
        '45',
        '5',
        '0',
        '5',
    ]
    dummio = DummyIO(to_write)
    get_new_config()


"""[
        'new_config',
        'Y',
        str(ALIGNMENTS_PATH),
        'fasta',
        '10, 150',
        '15, 25',
        '2',

        '2',
        'S1',
        '60, 61, 62, 80, 81, 82, 120',

        '2',
        'S2',
        '300, 301, 301',

        'ATCAGTCGATGCAGCTG',
        'ATCGATGCATGCTAGCT',

        'ACTGATGCATGCTAGTC',
        'ATCGATGCATGCTAGCT',

        '1',
        '3:4:5'

    ]"""