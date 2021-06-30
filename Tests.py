import pytest
from Bio.Seq import Seq

from hetero_spacer_generator import _gen_sequence_array, _select_and_set
from sequence_tools import *
from hetero_spacer_generator import *
import random

LARGE_SAMPLE_SIZE = 100
MED_SAMPLE_SIZE = 50
SMALL_SAMPLE_SIZE = 10

seqs = [
    Seq("ATCGATCG"),
    Seq("ATCG"),
    Seq("AAAAAAAA"),
    Seq("TTTTTTTT"),
    Seq(""),
    Seq("A")
]


def ensure_hetero_bases(bases: List[str]) -> bool:
    for i in range(len(bases)):
        for j in range(i + 1, len(bases)):
            if bases[i] == bases[j]:
                return False
    return True


def ensure_hetero_seq_arr(seq_arr: [List[List[str]]], num_hetero: int) -> bool:
    """Returns true iff the sequences contained in <seq_arr> are
    heterogeneous for the first <num_hetero> bases."""
    for column in range(min(num_hetero, len(seq_arr[0]) - 1)):
        bases = []
        for row in range(4):
            if seq_arr[row][column].isalpha():
                bases.append(seq_arr[row][column])
        if not ensure_hetero_bases(bases):
            return False
    return True


def ensure_valid_spacers(incomplete_primer: MBPrimerBuilder or Seq,
                         spacers: Tuple[int, int, int, int],
                         num_hetero) -> bool:
    """Ensures that the spacers produced by the given spacers will be valid. If
    """
    if type(incomplete_primer) == MBPrimerBuilder:
        seq_arr = _gen_sequence_array(incomplete_primer.get_binding_seq(),
                                      spacers)
    else:
        seq_arr = _gen_sequence_array(incomplete_primer,
                                      spacers)
    return ensure_hetero_seq_arr(seq_arr, num_hetero)


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


ADAPTER_LEN = 12
INDEX_LEN = 4
BINDING_MAX = 30
BINDING_MIN = 12


def gen_incomplete_primer(binding_len: int = 12):
    """Randomly generates an incomplete primer with a binding seq with len
    <binding_len> and blank heterogeneity regions. Sets the adapter and
    indexing regions to random sequences of default len."""
    primer = MBPrimerBuilder()
    primer.set_adapter_seq(gen_random_seq(ADAPTER_LEN))
    primer.set_binding_seq(gen_random_seq(binding_len))
    primer.set_index_seq(gen_random_seq(INDEX_LEN))
    return primer

def gen_rand_spacer_aligns(num_primers: int, max_val: int) \
        -> List[Tuple[int, int, int, int]]:
    """Generates a set of <num_primers> randomised primers with
    length <= <max_val>."""
    spacers = []
    for _ in range(num_primers):
        spacers.append((random.randrange(0, max_val),
                        random.randrange(0, max_val),
                        random.randrange(0, max_val),
                        random.randrange(0, max_val)))
    return spacers


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

        self.num_to_generate = LARGE_SAMPLE_SIZE

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
        self.incomplete_forward_primer = gen_incomplete_primer(forward_len)
        self.incomplete_reverse_primer = gen_incomplete_primer(reverse_len)

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


class TestSpacerAlignmentGen:

    def test_get_smallest_total_len_list(self):
        spacers = gen_rand_spacer_aligns(NUM_HETERO, NUM_HETERO)
        scores = get_smallest_total_len_list(spacers)

        for i in range(len(spacers)):
            assert sum(spacers[i]) == scores[i]

    def test_get_smallest_of_any_spacer_list(self) -> None:
        spacers = gen_rand_spacer_aligns(NUM_HETERO, NUM_HETERO)
        scores = get_smallest_of_any_spacer_list(spacers)

        for i in range(len(spacers)):
            assert max(spacers[i]) == scores[i]

    def test_sort_spacer_combos(self) -> None:
        spacers = gen_rand_spacer_aligns(NUM_HETERO, NUM_HETERO)
        sag = SpacerAlignmentGen(MAX_SPACER_LENGTH, NUM_HETERO)
        sag.sort_spacer_combos(spacers)
        assert sum(spacers[0]) < sum(spacers[-1]) and \
               max(spacers[0]) < max(spacers[-1])

    def test_get_all_spacer_combos_always_align(self):
        """When the heterogenity spacers are allowed to be as long as the
        heterogeneity region, there should always exist valid spacers."""
        sfm = SeqFixtureManager()
        sfm.set_primers_random()
        sfm.incomplete_forward_primer.set_binding_seq("AAAAAAAAAAAA")
        sag = SpacerAlignmentGen(12, 12)
        assert sag.get_all_spacer_combos(
            sfm.incomplete_forward_primer.get_binding_seq())

    def test_always_valid_align(self) -> None:
        for i in range(SMALL_SAMPLE_SIZE):
            for j in range(SMALL_SAMPLE_SIZE):
                sfm = SeqFixtureManager()
                sfm.for_num_hetero = i
                sfm.for_max_spacer_length = j
                sfm.set_primers_random()
                sfm.set_pot_spacers()
                for spacer in sfm.pot_forward_spacers:
                    assert ensure_valid_spacers(sfm.incomplete_forward_primer,
                                                spacer, i)
                for spacer in sfm.pot_reverse_spacers:
                    assert ensure_valid_spacers(sfm.incomplete_reverse_primer,
                                                spacer, i)

    def test_generates_valid_aligns(self):
        seq = Seq('GCCGGCATGGTCATGAAG')
        sag = SpacerAlignmentGen(NUM_HETERO, MAX_SPACER_LENGTH)
        primer_alignments = sag.get_all_spacer_combos(seq)
        for spacer in primer_alignments:
            assert ensure_valid_spacers(seq, spacer, NUM_HETERO)

    def test_is_invalid_seq(self):
        seq = Seq('GCCGGCATGGTCATGAAG')
        sag = SpacerAlignmentGen(NUM_HETERO, MAX_SPACER_LENGTH)
        for i in range(NUM_HETERO):
            for j in range(max(1, i), NUM_HETERO):
                assert not ensure_valid_spacers(seq, (0, 0, i, j), NUM_HETERO)


class TestRandomBaseSelection:
    """Test suite for methods that aid in the creation of random heterogeneity
    spacers"""

    random_tests_to_run = LARGE_SAMPLE_SIZE
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
            assert ensure_hetero_seq_arr(seq_arr, len(seq_arr[0]) - 1)

    def test_remove_high_dimer_complementarity_basic(self) -> None:
        """Tests whether _remove_high_dimer_complementarity produces a sample
        with the correct output size"""
        sfm = SeqFixtureManager()
        sfm.num_to_generate = MED_SAMPLE_SIZE
        sfm.do_all()
        rsg = RandomSpacerGen(12, 12, ConsolePresenter())
        rsg._remove_high_dimer_complementarity(sfm.forward_spacer_seqs,
                                               sfm.incomplete_forward_primer,
                                               SMALL_SAMPLE_SIZE)
        rsg._remove_high_dimer_complementarity(sfm.reverse_spacer_seqs,
                                               sfm.incomplete_reverse_primer,
                                               SMALL_SAMPLE_SIZE)
        assert len(sfm.forward_spacer_seqs) == \
               len(sfm.reverse_spacer_seqs) == SMALL_SAMPLE_SIZE

    def test_remove_high_dimer_complementarity_outprm_random(self) -> None:
        """Tests whether the _remove_high_dimer_complementarity produces primers
        with less binding than a random sampling."""
        sfm = SeqFixtureManager()
        sfm.num_to_generate = MED_SAMPLE_SIZE
        sfm.do_all()
        rsg = RandomSpacerGen(12, 12, ConsolePresenter())
        rand_for_spacers = []
        for i in range(SMALL_SAMPLE_SIZE):
            rand_for_spacers.append(random.choice(sfm.forward_spacer_seqs))
        rsg._remove_high_dimer_complementarity(sfm.forward_spacer_seqs,
                                               sfm.incomplete_forward_primer,
                                               SMALL_SAMPLE_SIZE)
        rand_scores = 0
        sorted_scores = 0
        for i in range(SMALL_SAMPLE_SIZE):
            sorted_scores += eval_total_complementarity(
                sfm.incomplete_forward_primer,
                sfm.forward_spacer_seqs[i])
            rand_scores += eval_total_complementarity(
                sfm.incomplete_forward_primer,
                rand_for_spacers[i])
        assert rand_scores >= sorted_scores

    def test_remove_high_consec_complementarity_basic(self) -> None:
        """Tests whether _remove_high_dimer_complementarity produces a sample with the
        correct output size"""
        sfm = SeqFixtureManager()
        sfm.num_to_generate = MED_SAMPLE_SIZE
        sfm.do_all()
        rsg = RandomSpacerGen(12, 12, ConsolePresenter())
        rsg._remove_high_consec_complementarity(sfm.forward_spacer_seqs,
                                                sfm.incomplete_forward_primer,
                                                SMALL_SAMPLE_SIZE)
        rsg._remove_high_consec_complementarity(sfm.reverse_spacer_seqs,
                                                sfm.incomplete_reverse_primer,
                                                SMALL_SAMPLE_SIZE)
        assert len(sfm.forward_spacer_seqs) == \
               len(sfm.reverse_spacer_seqs) == SMALL_SAMPLE_SIZE

    def test_filter_spacer_sets_basic(self) -> None:
        sfm = SeqFixtureManager()
        sfm.num_to_generate = MED_SAMPLE_SIZE
        sfm.do_all()
        rsg = RandomSpacerGen(12, 12, ConsolePresenter())
        rsg._filter_spacer_sets(sfm.incomplete_forward_primer,
                                sfm.incomplete_reverse_primer,
                                sfm.forward_spacer_seqs,
                                sfm.reverse_spacer_seqs,
                                SMALL_SAMPLE_SIZE / MED_SAMPLE_SIZE * 100)
        assert len(sfm.forward_spacer_seqs) == \
               len(sfm.reverse_spacer_seqs) == SMALL_SAMPLE_SIZE

    def test_cross_comapre_basic(self):
        """Tests that cross compare prodcues the expected output for a basic
        input."""
        sfm = SeqFixtureManager()
        sfm.do_all()
        rsg = RandomSpacerGen(12, 12)
        numsets = 3
        primer_sets = rsg._cross_compare(sfm.incomplete_forward_primer,
                                         sfm.incomplete_reverse_primer,
                                         sfm.forward_spacer_seqs[
                                         0:SMALL_SAMPLE_SIZE],
                                         sfm.reverse_spacer_seqs[
                                         0:SMALL_SAMPLE_SIZE],
                                         numsets)
        assert len(primer_sets) == numsets

    def test_full_primer_creation_process(self):
        for_adapter_seq = 'ACACTCTTTCCCTACACGACGCTCTTCCGATCT'
        rev_adapter_seq = 'GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT'
        for_indexing_seq = 'ATCG'
        rev_indexing_seq = 'GCTA'
        for_binding_seq = 'TGTCATCTCCTTCTGCTACGG'
        rev_binding_seq = 'GCCGGCATGGTCATGAAG'

        hetero_size = 12
        spacer_size = 12

        incomp_forward_primer = MBPrimerBuilder()
        incomp_forward_primer.set_binding_seq(for_binding_seq)
        incomp_forward_primer.set_index_seq(for_indexing_seq)
        incomp_forward_primer.set_adapter_seq(for_adapter_seq)
        incomp_reverse_primer = MBPrimerBuilder()
        incomp_reverse_primer.set_binding_seq(rev_binding_seq)
        incomp_reverse_primer.set_index_seq(rev_indexing_seq)
        incomp_reverse_primer.set_adapter_seq(rev_adapter_seq)

        hg = HeteroGen(hetero_size, spacer_size)
        hg.set_rigour(-5)
        for_spacers = hg.get_all_spacer_combos(incomp_forward_primer
                                               .get_binding_seq())
        rev_spacers = hg.get_all_spacer_combos(incomp_reverse_primer
                                               .get_binding_seq())

        primers = hg.get_hetero_seqs(incomp_forward_primer,
                                     incomp_reverse_primer,
                                     for_spacers[1],
                                     rev_spacers[1],
                                     5)
        assert len(primers) == 5


"""
    def test_remove_high_dimerisation(self) -> None:

    def test_cross_comapre(self) -> None:
    """

if __name__ == '__main__':
    pytest.main(['Tests.py'])
