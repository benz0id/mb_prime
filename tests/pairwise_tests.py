import random

import pytest
from Bio.Seq import Seq

import fixtures_and_helpers as fah
from hetero_spacer_generator.sequence_tools import SeqAnalyzer


class TestPairwiseSorter:
    sfm = fah.SeqFixtureManager()
    sfm.num_to_generate = fah.MED_SAMPLE_SIZE
    sfm.do_all()

    pass



class TestSeqAnalyser:

    def test_com_seqs_any_overlap_prod_seqs_same_len(self) -> None:

        def equal_length(seq1: str, seq2: str) -> int:
            assert len(seq1) == len(seq2)
            return 1

        seqa = SeqAnalyzer()

        for i in range(fah.SMALL_SAMPLE_SIZE):
            seq1 = fah.gen_random_seq(random.randrange(0, fah.SMALL_SEQ_LEN))
            seq2 = fah.gen_random_seq(random.randrange(0, fah.SMALL_SEQ_LEN))
            seqa.comp_seqs_any_overlap(seq1, seq2, equal_length)

    def test_consec_complementarity(self) -> None:

        seqa = SeqAnalyzer()
        seqa.expect_degeneracy(False)

        seq1 = Seq('ACAAGTG CTAGT AGCTGATCG')  # 5' - 3'
        seq2 = Seq('TTGATGC GATCA ATGCTGTAG')  # 3' - 5'
        # Matching          ^^^^^

        assert seqa.comp_seqs_any_overlap(seq1, seq2,
                                          seqa.get_consec_complementarity) == 5
"""
if __name__ == '__main__':
    pytest.main(['tests/pairwise_tests', '-v'])"""
