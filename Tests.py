import pytest
from Bio.Seq import Seq

from hetero_spacer_generator import _gen_sequence_array, _select_and_set
from sequence_tools import *
from hetero_spacer_generator import *
import random

seqs = [
    Seq("ATCGATCG"),
    Seq("ATCG"),
    Seq("AAAAAAAA"),
    Seq("TTTTTTTT"),
    Seq(""),
    Seq("A")
]


def ensure_hetero_bases(bases: List[str]) -> bool:
    for i in range(4):
        for j in range(i + 1, 4):
            if bases[i] == bases[j]:
                return False
    return True


def ensure_hetero_seq_arr(seq_arr: [List[List[str]]], num_hetero: int) -> bool:
    """Returns true iff the sequences contained in <seq_arr> are
    heterogeneous for the first <num_hetero> bases."""
    for column in range(num_hetero):
        bases = []
        for row in range(4):
            bases.append(seq_arr[row][column])
        if not ensure_hetero_bases(bases):
            return False
    return True


def gen_random_seq(length: int) -> Seq:
    """Generates a random Seq of <length>"""
    bases = ['A', 'T', 'C', 'G']
    seq_str = ''
    for _ in range(length):
        seq_str += random.choice(bases)

    return Seq(seq_str)


def gen_random_spacers(incomplete_forward_primer: MBPrimerBuilder,
                       incomplete_reverse_primer: MBPrimerBuilder,
                       forward_spacer_seqs: List[Tuple[Seq]],
                       reverse_spacer_seqs: List[Tuple[Seq]],
                       forward_spacer: Tuple[int, int, int, int],
                       reverse_spacer: Tuple[int, int, int, int],
                       num_to_generate: float) -> None:
    """Empties <forward_spacer_seqs> and <reverse_spacer_seqs>, refilling them
    with <num_to_generate> randomly generated sets of heterogeneity spacers.
     Generated spacers with be based upon <incomplete_forward_primer> and
     <incomplete_reverse_primer>."""
    forward_spacer_seqs.clear()
    reverse_spacer_seqs.clear()
    generic_rsg = RandomSpacerGen(12, 12, ConsolePresenter())
    generic_rsg._random_per_align = num_to_generate
    forward_spacers = generic_rsg._gen_hetero_set(incomplete_forward_primer,
                                                  forward_spacer)
    forward_spacer_seqs.extend(forward_spacers)
    reverse_spacers = generic_rsg._gen_hetero_set(incomplete_reverse_primer,
                                                  reverse_spacer)
    reverse_spacer_seqs.extend(reverse_spacers)


class TestSequenceTools:
    """Test suite for sequence tools"""

    def test_get_max_complementarity(self):
        seq = seqs[1]
        assert get_max_complementarity(seq, [seqs[0]]) == 2
        assert get_max_complementarity(seq, [seqs[1]]) == 0
        assert get_max_complementarity(seq, [seqs[2]]) == 1
        assert get_max_complementarity(seq, [seqs[3]]) == 1
        assert get_max_complementarity(seq, [seqs[4]]) == 0
        assert get_max_complementarity(seq, [seqs[5]]) == 1

    def test_get_site_complementarity(self):
        assert get_site_complementarity(seqs[4], seqs[1], 1) == 0
        assert get_site_complementarity(seqs[1], seqs[0], 0) == 0
        assert get_site_complementarity(seqs[1], seqs[0], 1) == 2


def test_co_sort() -> None:
    ls = [1, 2, 3, 4, 5, 6, 7]
    lst = [2, 7, 1, 3, 5, 6, 4]
    lst1 = [2, 7, 1, 3, 5, 6, 4]
    lstc = lst[:]
    lstc1 = lst1[:]

    co_sort(lstc, lstc1)
    assert lstc == ls and lstc1 == ls


ADAPTER_LEN = 12
INDEX_LEN = 4
BINDING_MAX = 30
BINDING_MIN = 12


def gen_incomplete_primers(binding_len: int = 12):
    """Randomly generates an incomplete primer with a binding seq with len
    <binding_len> and blank heterogeneity regions. Sets the adapter and
    indexing regions to random sequences of default len."""
    primer = MBPrimerBuilder()
    primer.set_adapter_seq(gen_random_seq(ADAPTER_LEN))
    primer.set_binding_seq(gen_random_seq(binding_len))
    primer.set_index_seq(gen_random_seq(INDEX_LEN))
    return primer


class SeqFixtureManager:
    """Manages various datastructures required for testing."""
    incomplete_forward_primer: MBPrimerBuilder
    incomplete_reverse_primer: MBPrimerBuilder

    for_num_hetero: int
    rev_num_hetero: int

    for_max_spacer_length: int
    rev_max_spacer_length: int

    pot_forward_spacers: List[Tuple[int, int, int, int]]
    pot_reverse_spacers: List[Tuple[int, int, int, int]]

    forward_spacer: Tuple[int, int, int, int]
    reverse_spacer: Tuple[int, int, int, int]

    forward_spacer_seqs = []
    reverse_spacer_seqs = []

    num_to_generate: int

    def __init__(self) -> None:
        self.rsg4 = RandomSpacerGen(4, 4, ConsolePresenter())
        self.rsg12 = RandomSpacerGen(12, 12, ConsolePresenter())
        self.rsg_12_high_rigour = RandomSpacerGen(12, 12, ConsolePresenter(), 1)

        self.incomplete_forward_primer = MBPrimerBuilder()
        self.incomplete_reverse_primer = MBPrimerBuilder()

        self.for_num_hetero = 12
        self.rev_num_hetero = 12

        self.for_max_spacer_length = 12
        self.rev_max_spacer_length = 12

        self.pot_forward_spacers = []
        self.pot_reverse_spacers = []

        self.forward_spacer = (1, 1, 1, 1)
        self.reverse_spacer = (1, 1, 1, 1)

        self.forward_spacer_seqs = []
        self.reverse_spacer_seqs = []

        self.num_to_generate = 1000

    def do_all(self):
        """Constructs all contained attributes."""
        self.set_primers_random()
        self.set_pot_spacers()
        self.set_single_spacers_random()
        self.set_spacers_seqs_random()

    def set_primers_random(self):
        """Sets the <incomplete_forward_primer> and <incomplete_reverse_primer>
        to random primers with region length generated as specified by defaults.
        """
        forward_len = random.randrange(BINDING_MIN, BINDING_MAX)
        reverse_len = random.randrange(BINDING_MIN, BINDING_MAX)
        self.incomplete_forward_primer = gen_incomplete_primers(forward_len)
        self.incomplete_reverse_primer = gen_incomplete_primers(reverse_len)

    def set_pot_spacers(self):
        """Sets <pot_forward_spacer> and <pot_reverse_spacer> to a set of
        heterogeneity spacer alignments."""
        sag = SpacerAlignmentGen(self.for_max_spacer_length,
                                 self.for_num_hetero, ConsolePresenter())
        self.pot_forward_spacers = sag.get_all_spacer_combos(
            self.incomplete_forward_primer.get_binding_seq())

        sag.set_params(self.rev_max_spacer_length, self.rev_num_hetero)
        self.pot_reverse_spacers = sag.get_all_spacer_combos(
            self.incomplete_reverse_primer.get_binding_seq())

    def set_single_spacers_random(self):
        """Sets <self.forward_spacer> and <self.reverse_spacer> to a random
        spacer from the forward and reverse spacers."""
        self.forward_spacer = random.choice(self.pot_forward_spacers)
        self.reverse_spacer = random.choice(self.pot_reverse_spacers)

    def set_spacers_seqs_random(self):
        """Generates random forward and reverse spacers seqs using the
        attributes above."""
        gen_random_spacers(self.incomplete_forward_primer,
                           self.incomplete_reverse_primer,
                           self.forward_spacer_seqs,
                           self.reverse_spacer_seqs,
                           self.forward_spacer,
                           self.reverse_spacer,
                           self.num_to_generate)


class TestRandomBaseSelection:
    """Test suite for methods that aid in the creation of random heterogeneity
    spacers"""

    random_tests_to_run = 500
    hg = HeteroGen(num_hetero=3)
    hg1 = HeteroGen()
    simple_seq_arr = [['A', 'T', 'C', 'G'],
                      ['', 'A', 'T', 'C'],
                      ['', '', 'A', 'T'],
                      ['', '', '', 'A']]

    def test_gen_sequence_array(self) -> None:

        assert _gen_sequence_array(Seq('ATCG'), (0, 1, 2, 3)) \
               == self.simple_seq_arr

        for _ in range(self.random_tests_to_run):
            seq = gen_random_seq(12)
            spacers = self.hg1._alignment_gen.get_all_spacer_combos(seq)
            for spacer in spacers:
                _gen_sequence_array(seq, spacer)

    def test_get_vacant_bases(self) -> None:
        assert self.hg._primer_gen._get_vacant_bases((0, 1, 2, 3)) == [
            [1, 2, 3], [2, 3],
            [3], []]

        for _ in range(self.random_tests_to_run):
            seq = gen_random_seq(12)
            spacers = self.hg1._alignment_gen.get_all_spacer_combos(seq)
            for spacer in spacers:
                self.hg1._primer_gen._get_vacant_bases(spacer)

    def test_get_potential_bases(self):

        for i in range(4):
            pbs = self.hg._primer_gen._get_potential_bases(self.simple_seq_arr,
                                                           i)
            for j in range(4):
                assert self.simple_seq_arr[j][i] not in pbs

    def test_select_and_set_simple(self):

        potential_bases = ['T', 'C', 'G']
        unfilled_bases = [[1, 2, 3], []]
        column = 0
        seq_arr = self.simple_seq_arr.copy()
        for i in range(3):
            _select_and_set(potential_bases, unfilled_bases,
                            column, seq_arr)
            assert len(potential_bases) == 2 - i
            assert len(unfilled_bases[0]) == 2 - i
            for row in range(1, 4):
                assert seq_arr[row][column] != 'A'

    def test_gen_heterogeneity_spacer_rand(self):
        seq = gen_random_seq(12)
        spacer_combos = self.hg1._alignment_gen.get_all_spacer_combos(seq)
        if not spacer_combos:
            self.test_gen_heterogeneity_spacer_rand()
            return
        for spacer_lengths in spacer_combos:
            spacers = self.hg1._primer_gen._gen_heterogeneity_spacers_rand(seq,
                                                                spacer_lengths)
            seq_arr = _gen_sequence_array(seq, spacer_lengths)
            for i in range(4):
                for j in range(len(spacers[i])):
                    seq_arr[i][j] = spacers[i][j]
            assert ensure_hetero_seq_arr(seq_arr, 12)

    def test_filter_spacer_sets_basic(self) -> None:
        sfm = SeqFixtureManager()
        sfm.num_to_generate = 1000
        sfm.do_all()
        rsg = RandomSpacerGen(12, 12, ConsolePresenter())
        rsg._filter_spacer_sets(sfm.incomplete_forward_primer,
                                sfm.incomplete_reverse_primer,
                                sfm.forward_spacer_seqs,
                                sfm.reverse_spacer_seqs, 3)
        assert len(sfm.forward_spacer_seqs) == \
               len(sfm.reverse_spacer_seqs) == 30


"""
    def test_remove_high_dimerisation(self) -> None:

    def test_cross_comapre(self) -> None:
    """

if __name__ == '__main__':
    pytest.main(['Tests.py'])
