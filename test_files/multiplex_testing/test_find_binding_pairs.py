import pytest as pt

from config_handling.formatting import *
from seq_alignment_analyser.find_binding_pairs import *
from test_files.fixtures_and_helpers import ALIGNMENTS_PATH, configure_log_out

configure_log_out('FindBindingPairs')

def test_initialisation() -> None:
    target_sites = [
        TargetRegionInfo('test_target_1', 'alex_test_align.fas', [250, 260])
    ]

    adapters = [
        AdapterPair('TAGCTAGCTAGCTAGCTAGGG', 'ACGATGCATGTTGCTAGTGT')
    ]

    primer_params = PrimerParams(10, 10, InclRange(100, 200), InclRange(15, 25),
                                 20)

    alignments_path = ALIGNMENTS_PATH

    targ_mt = 65

    max_mt_deviance = 5

    fbp = FindBindingPairs(target_sites, adapters, primer_params,
                           alignments_path, targ_mt, max_mt_deviance)

    lol = 'dummy'


@pt.fixture
def single_target_fbp() -> FindBindingPairs:
    target_sites = [
        TargetRegionInfo('test_target_1', 'alex_test_align.fas', [250, 260])
    ]

    adapters = [
        AdapterPair('TAGCTAGCTAGCTAGCTAGGG', 'ACGATGCATGTTGCTAGTGT')
    ]

    primer_params = PrimerParams(10, 10, InclRange(100, 200), InclRange(15, 25),
                                 20)

    alignments_path = ALIGNMENTS_PATH

    targ_mt = 65

    max_mt_deviance = 5

    return FindBindingPairs(target_sites, adapters, primer_params,
                           alignments_path, targ_mt, max_mt_deviance)


def test_single_target_fbp(single_target_fbp: FindBindingPairs) -> None:
    single_target_fbp.get_n_most_conserved_valid(
        1000, single_target_fbp._iterator_manager.__next__()
    )