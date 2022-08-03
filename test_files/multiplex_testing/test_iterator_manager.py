import collections
from typing import Any
import logging
import pytest as pt
from seq_alignment_analyser.iterator_manager import *
import datetime
import pytest as pt
from test_files.fixtures_and_helpers import get_msa, ALIGNMENTS_PATH, \
MULTIPLEX_TESTING_PATH

handler = logging.FileHandler(MULTIPLEX_TESTING_PATH /
                              'IteratorManagerTest.log')
root_log = logging.getLogger('root')
root_log.addHandler(handler)
root_log.setLevel(10)
root_log.info('\n   ====   Test Session Begins: ' + str(datetime.datetime.now())
              + '   ====')




BIMinit = collections.namedtuple(
    'BIMinit',
    [
        'msa_to_targets',  'primer_pool', 'primer_primer_distance',
        'primer_target_distance', 'target_region_len', 'binding_region_len',
        'ideal_binding_size'
    ])


@pt.fixture
def basic_config() -> BIMinit:
    alex_msa = get_msa('alex_test_align')
    primer_primer_distance = 0
    primer_target_distance = 0
    target_region_len = InclRange(5, 5)
    binding_region_len = InclRange(1, 1)
    ideal_binding_size = 1
    target_1 = TargetRegionInfo('target_1', 'alex_test_align.fas',
                                [200])
    return BIMinit({alex_msa: [target_1]}, [], primer_primer_distance,
                   primer_target_distance, target_region_len,
                   binding_region_len, ideal_binding_size)


def test_BIM_functionality(basic_config) -> None:
    """Tests the basic functionality of the BindingIteratorManager"""
    bim = BindingIteratorManager(*basic_config)
    iterator = bim.iterator_queue.pop()
    num_pos = iterator.get_num_pos_primers()
    t_num_pos = 0
    for _ in iterator:
        t_num_pos += 1

    assert num_pos == t_num_pos == 5


@pt.fixture
def complex_config() -> BIMinit:
    alex_msa = get_msa('alex_test_align')
    primer_primer_distance = 10
    primer_target_distance = 20
    target_region_len = InclRange(60, 100)
    binding_region_len = InclRange(16, 24)
    ideal_binding_size = 1
    target_1 = TargetRegionInfo('target_1', 'alex_test_align.fas',
                                [200, 210])
    return BIMinit({alex_msa: [target_1]}, [], primer_primer_distance,
                   primer_target_distance, target_region_len,
                   binding_region_len, ideal_binding_size)


def test_complex_BIM_functionality(complex_config) -> None:
    """Tests the basic functionality of the BindingIteratorManager"""
    bim = BindingIteratorManager(*complex_config)
    iterator = bim.iterator_queue.pop()
    num_pos = iterator.get_num_pos_primers()
    t_num_pos = 0
    for _ in iterator:
        t_num_pos += 1

    assert num_pos == t_num_pos > 30


@pt.fixture
def basic_overlap_config() -> BIMinit:
    alex_msa = get_msa('alex_test_align')
    primer_primer_distance = 0
    primer_target_distance = 0
    target_region_len = InclRange(2, 6)
    binding_region_len = InclRange(1, 1)
    ideal_binding_size = 1

    target_1 = TargetRegionInfo('target_1', 'alex_test_align.fas',
                                [20, 21])
    target_2 = TargetRegionInfo('target_2', 'alex_test_align.fas',
                                [24, 25])
    return BIMinit({alex_msa: [target_1, target_2]}, [], primer_primer_distance,
                   primer_target_distance, target_region_len,
                   binding_region_len, ideal_binding_size)


def test_basic_overlap_BIM_functionality(basic_overlap_config) -> None:
    """Tests the basic functionality of the BindingIteratorManager"""
    bim = BindingIteratorManager(*basic_overlap_config)
    iterator1 = bim.iterator_queue.pop()
    iterator2 = bim.iterator_queue.pop()

    assert iterator1.get_num_pos_primers() == \
        iterator2.get_num_pos_primers() == 5


@pt.fixture
def complex_overlap_config() -> BIMinit:
    alex_msa = get_msa('alex_test_align')
    primer_primer_distance = 10
    primer_target_distance = 10
    target_region_len = InclRange(0, 100)
    binding_region_len = InclRange(20, 20)
    ideal_binding_size = 20

    target_1 = TargetRegionInfo('target_1', 'alex_test_align.fas',
                                [140])
    target_2 = TargetRegionInfo('target_2', 'alex_test_align.fas',
                                [211])
    target_3 = TargetRegionInfo('target_3', 'alex_test_align.fas',
                                [282])

    print(alex_msa.__len__())
    return BIMinit({alex_msa: [target_1, target_2, target_3]}, [], 
                   primer_primer_distance, primer_target_distance,
                   target_region_len, binding_region_len, ideal_binding_size)


def test_complex_overlap_BIM_functionality(complex_overlap_config) -> None:
    """Tests the basic functionality of the BindingIteratorManager"""

    bim = BindingIteratorManager(*complex_overlap_config)
    iterator1 = bim.iterator_queue.popleft()
    iterator2 = bim.iterator_queue.popleft()
    iterator3 = bim.iterator_queue.popleft()

    iterator1_size = iterator1.get_num_pos_primers()
    iterator2_size = iterator2.get_num_pos_primers()
    iterator3_size = iterator3.get_num_pos_primers()

    assert iterator1_size <= iterator2_size <= iterator3_size








