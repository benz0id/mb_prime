import pytest
from Bio.Seq import Seq
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
    """Generates a ranod Seq of <length>"""
    bases = ['A', 'T', 'C', 'G']
    seq_str = ''
    for _ in range(length):
        seq_str += random.choice(bases)

    return Seq(seq_str)


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

        assert self.hg._gen_sequence_array(Seq('ATCG'), (0, 1, 2, 3)) \
               == self.simple_seq_arr

        for _ in range(self.random_tests_to_run):
            seq = gen_random_seq(12)
            spacers = self.hg1.get_all_spacer_combos(seq)
            for spacer in spacers:
                self.hg1._gen_sequence_array(seq, spacer)

    def test_get_vacant_bases(self) -> None:
        assert self.hg.get_vacant_bases((0, 1, 2, 3)) == [[1, 2, 3], [2, 3],
                                                          [3], []]

        for _ in range(self.random_tests_to_run):
            seq = gen_random_seq(12)
            spacers = self.hg1.get_all_spacer_combos(seq)
            for spacer in spacers:
                self.hg1.get_vacant_bases(spacer)

    def test_get_potential_bases(self):

        for i in range(4):
            pbs = self.hg.get_potential_bases(self.simple_seq_arr, i)
            for j in range(4):
                assert self.simple_seq_arr[j][i] not in pbs

    def test_select_and_set_simple(self):

        potential_bases = ['T', 'C', 'G']
        unfilled_bases = [[1, 2, 3], []]
        column = 0
        seq_arr = self.simple_seq_arr.copy()
        for i in range(3):
            self.hg.select_and_set(potential_bases, unfilled_bases,
                                   column, seq_arr)
            assert len(potential_bases) == 2 - i
            assert len(unfilled_bases[0]) == 2 - i
            for row in range(1, 4):
                assert seq_arr[row][column] != 'A'

    def test_gen_heterogeneity_spacer_rand(self):
        seq = gen_random_seq(12)
        spacer_combos = self.hg1.get_all_spacer_combos(seq)
        if not spacer_combos:
            self.test_gen_heterogeneity_spacer_rand()
            return
        for spacer_lengths in spacer_combos:
            spacers = self.hg1.gen_heterogeneity_spacers_rand(seq,
                                                              spacer_lengths)
            seq_arr = self.hg1._gen_sequence_array(seq, spacer_lengths)
            for i in range(4):
                for j in range(len(spacers[i])):
                    seq_arr[i][j] = spacers[i][j]
            self.hg1.visualise_seq_arr([seq_arr])
            assert ensure_hetero_seq_arr(seq_arr, 12)


if __name__ == '__main__':
    pytest.main(['Tests.py'])
