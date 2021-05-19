from Bio.Seq import Seq
from typing import Tuple
from primer_tools import MB_Primer


def generate_hetero_seq(primer: Seq, vary_len: int,
                        num_sets: int) -> Tuple[Seq]:
    """Generates heterogeneity sequences given a <primer> sequence. Ensures
    nucleotide diversity among the first <vary_len> bases in the sequence.
    Returns a list of <num_sets> tuples, each containing 4 tandem heterogeneity
    sequences.
    Precondition: len(primer) > vary_len."""

    seq1

def align_hetero(seq: Seq, num_hetero: int) -> Tuple[int]:
    """Given a <seq>, will find 4 alignments of that sequence (produced by
    shifting it to the right) that ensures nucleotide diversity across the first
    <num_hetero> bases.

    >>> seq = Seq('ATCGAA')
    >>> align_hetero(seq, 4)
    (0, 1, 2, 3)

    Corresponding to the alignment:
    ATCGAA
    +ATCGAA
    ++ATCGAA
    +++ATCGAA

    >>> seq = Seq('ATCGAA')
    >>> align_hetero(seq, 5)




    """
