import sys
import random
import pytest as pt

from config_handling.formatting import *
from seq_alignment_analyser.find_binding_pairs import *
from test_files.fixtures_and_helpers import ALIGNMENTS_PATH, configure_log_out, \
    get_msa

configure_log_out('FindBindingPairs')

@pt.fixture
def single_target_fbp() -> FindBindingPairs:
    target_sites = [
        TargetRegionInfo('test_target_1', 'alex_test_align.fas', [250, 260])
    ]

    adapters = [
        AdapterPair('TAGCTAGCTAGCTAGCTAGGG', 'ACGATGCATGTTGCTAGTGT')
    ]

    primer_params = PrimerParams(10, 10, InclRange(70, 200), InclRange(15, 25),
                                 20, 126)

    alignments_path = ALIGNMENTS_PATH

    targ_mt = 60

    max_mt_deviance = 5

    fbp = FindBindingPairs(target_sites, adapters, primer_params,
                           alignments_path, targ_mt, max_mt_deviance)
    fbp._cur_msa = get_msa('alex_test_align')
    fbp._cur_target = 'test_target_1'
    return fbp


def test_single_target_fbp(single_target_fbp: FindBindingPairs) -> None:
    single_target_fbp.get_n_most_conserved_valid(
        1000, single_target_fbp._iterator_manager.__next__())

def test_single_target_fbp_half(single_target_fbp: FindBindingPairs) -> None:
    single_target_fbp.get_best_binding_pair(
        single_target_fbp._iterator_manager.__next__())

@pt.fixture
def multi_target_fbp_no_conflict_open() -> FindBindingPairs:
    target_sites = [
        TargetRegionInfo('test_target_1', 'alex_test_align.fas', [250, 260]),
        TargetRegionInfo('test_target_2', 'alex_test_align.fas', [100, 110])
    ]

    adapters = [
        AdapterPair('TAGCTAGCTAGCTAGCTAGGG', 'ACGATGCATGTTGCTAGTGT'),
        AdapterPair('TAGCTAGCTAGCTAGCTAGGG', 'ACGATGCATGTTGCTAGTGT')
    ]

    primer_params = PrimerParams(10, 10, InclRange(70, 200), InclRange(15, 25),
                                 20, 126)

    alignments_path = ALIGNMENTS_PATH

    targ_mt = 60

    max_mt_deviance = 5

    return FindBindingPairs(target_sites, adapters, primer_params,
                           alignments_path, targ_mt, max_mt_deviance, mode=OPEN)

def test_multi_target_fbp_no_conflict(multi_target_fbp_no_conflict_open: FindBindingPairs) -> None:
    multi_target_fbp_no_conflict_open.get_best_binding_pairs()

@pt.fixture
def multi_target_fbp_no_conflict_restricted() -> FindBindingPairs:
    target_sites = [
        TargetRegionInfo('test_target_1', 'alex_test_align.fas', [250, 260]),
        TargetRegionInfo('test_target_2', 'alex_test_align.fas', [100, 110])
    ]

    adapters = [
        AdapterPair('TAGCTAGCTAGCTAGCTAGGG', 'ACGATGCATGTTGCTAGTGT'),
        AdapterPair('TAGCTAGCTAGCTAGCTAGGG', 'ACGATGCATGTTGCTAGTGT')
    ]

    primer_params = PrimerParams(10, 10, InclRange(70, 200), InclRange(15, 25),
                                 20, 126)

    alignments_path = ALIGNMENTS_PATH

    targ_mt = 60

    max_mt_deviance = 5

    return FindBindingPairs(target_sites, adapters, primer_params,
                           alignments_path, targ_mt, max_mt_deviance, mode=RESTRICTED)

def test_multi_target_fbp_no_conflict_restricted(
        multi_target_fbp_no_conflict_restricted: FindBindingPairs) -> None:
    multi_target_fbp_no_conflict_restricted.get_best_binding_pairs()


@pt.fixture
def multi_target_fbp_no_conflict_complex() -> FindBindingPairs:
    target_sites = [
        TargetRegionInfo('test_target_1', 'alex_test_align.fas', [100, 110]),
        TargetRegionInfo('test_target_2', 'alex_test_align.fas', [200, 210]),
        TargetRegionInfo('test_target_3', 'alex_test_align.fas', [300, 310]),
        TargetRegionInfo('test_target_4', 'alex_test_align.fas', [400, 410])
    ]

    adapters = [
        AdapterPair('TAGCTAGCTAGCTAGCTAGGG', 'ACGATGCATGTTGCTAGTGT'),
        AdapterPair('AGCTAGCTAGCTAGTTACGGG', 'CAAACTAGCTAGTGTGGGGA'),
        AdapterPair('CGTAGCTACGTAGCTAGCTGA', 'AGTCATCGATGCTAGCTAGA'),
        AdapterPair('TAGCTAGCTAGCTGASTCGAT', 'ACGTGATGCATGTACTGATC')
    ]

    primer_params = PrimerParams(10, 0, InclRange(40, 150), InclRange(15, 25),
                                 20, 126)

    alignments_path = ALIGNMENTS_PATH

    targ_mt = 60

    max_mt_deviance = 5

    return FindBindingPairs(target_sites, adapters, primer_params,
                           alignments_path, targ_mt, max_mt_deviance)

def test_multi_target_fbp_no_conflict_complex(multi_target_fbp_no_conflict_complex) -> None:
    multi_target_fbp_no_conflict_complex.get_best_binding_pairs()


@pt.fixture
def multi_target_fbp_no_conflict_complex_open() -> FindBindingPairs:
    target_sites = [
        TargetRegionInfo('test_target_1', 'alex_test_align.fas', [100, 110]),
        TargetRegionInfo('test_target_2', 'alex_test_align.fas', [200, 210]),
        TargetRegionInfo('test_target_3', 'alex_test_align.fas', [300, 310]),
        TargetRegionInfo('test_target_4', 'alex_test_align.fas', [400, 410])
    ]

    adapters = [
        AdapterPair('TAGCTAGCTAGCTAGCTAGGG', 'ACGATGCATGTTGCTAGTGT'),
        AdapterPair('AGCTAGCTAGCTAGTTACGGG', 'CAAACTAGCTAGTGTGGGGA'),
        AdapterPair('CGTAGCTACGTAGCTAGCTGA', 'AGTCATCGATGCTAGCTAGA'),
        AdapterPair('TAGCTAGCTAGCTGASTCGAT', 'ACGTGATGCATGTACTGATC')
    ]

    primer_params = PrimerParams(10, 10, InclRange(20, 150), InclRange(15, 25),
                                 20, 126)

    alignments_path = ALIGNMENTS_PATH

    targ_mt = 60

    max_mt_deviance = 5

    return FindBindingPairs(target_sites, adapters, primer_params,
                           alignments_path, targ_mt, max_mt_deviance)

def test_multi_target_fbp_no_conflict_complex_open(multi_target_fbp_no_conflict_complex_open) -> None:
    multi_target_fbp_no_conflict_complex_open.get_best_binding_pairs()

rand_seq = lambda n: ''.join([random.choice('ATCG-') for _ in range(n)])

def test_bulk() -> None:
    """msa_len = 2000
    num_seqs = 5

    seqs = [rand_seq(msa_len) for _ in range(num_seqs)]

    msa = MSA(ALIGNMENTS_PATH / 'alex_test_align.fas')
    msa._seqs = seqs
    msa._parse_consensus_attributes()



    while next_targ"""


