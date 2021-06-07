import pytest
from Bio.Seq import Seq
from sequence_tools import *

seqs = [
    Seq("ATCGATCG"),
    Seq("ATCG"),
    Seq("AAAAAAAA"),
    Seq("TTTTTTTT"),
    Seq(""),
    Seq("A")
]


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


if __name__ == '__main__':
    pytest.main(['Tests.py'])
