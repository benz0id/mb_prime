from Bio.Seq import Seq
from typing import List, Dict, Iterable, Callable


def iterate_seqs(oligo: Seq, seqs: Iterable[Seq],
                 comp_method: Callable[[Seq, Seq, int], int]) -> int:
    """Compares all completely overlapping combinations of <oligo> over each of
     <seq> using <comp_method>, returns the greatest complementarity found."""
    max_comp = 0
    oligo_len = len(oligo)
    for seq in seqs:
        # Each possible alignment along seq
        seq_len = len(seq)
        if seq_len >= oligo_len:  # Iterate over alignments of oligo over seq.
            larger_seq = seq
            smaller_seq = oligo
            num_alignments = len(seq) - oligo_len
        else:  # Iterate over alignments of seq over oligo.
            larger_seq = oligo
            smaller_seq = seq
            num_alignments = oligo_len - len(seq)

        for offset in range(0, num_alignments):
            local_max_comp = comp_method(smaller_seq,
                                                      larger_seq, offset)
            if local_max_comp > max_comp:
                max_comp = local_max_comp
    return max_comp


def get_max_complementarity(oligo: Seq, seqs: Iterable[Seq]) -> int:
    """Returns the highest number of complimentary bases between <oligo> at any
    possible alignment on any of <seqs>, while <oligo> and <seqs> are still
    fully overlapped.

    >>> get_max_complementarity(Seq("ATCGA"), [Seq("ATTCAGACC")])
     4"""
    return iterate_seqs(oligo, seqs, get_site_complementarity)


def get_site_complementarity(oligo: Seq, seq: Seq, offset: int) -> int:
    """Gets the number of complementary bases of <oligo> at the site along
    <seq> specified by  the <offset>.

    Precondition:
    offset + len(oligo) <= len(seq)
    >>> get_site_complementarity(Seq("ATCG"), Seq("TTAG CC"), 0)
     1
    >>> get_site_complementarity(Seq("ATCG"), Seq("T TAGC C"), 1)
     4
     """
    num_complement = 0
    oligo_compliment = oligo.complement()

    # At each base, check for equality between the oligo's compliment and seq
    for i in range(len(oligo)):
        if oligo_compliment[i] == seq[i + offset]:
            num_complement += 1
    return num_complement


def get_max_complementarity_consec(oligo: Seq, seqs: Iterable[Seq]) -> int:
    """Returns the highest number of consicutive complimentary bases between
    <oligo> at any possible alignment on any of <seqs>, while <oligo> and <seqs>
    are still fully overlapped.

    >>> get_max_complementarity(Seq("ATCGA"), [Seq("ATTCAGACC")])
     4"""
    return iterate_seqs(oligo, seqs, get_site_complementarity_consec)


def get_site_complementarity_consec(oligo: Seq, seq: Seq, offset: int) -> int:
    """Gets the number of complementary bases of <oligo> at the site along
    <seq> specified by  the <offset>.

    Precondition:
    offset + len(oligo) <= len(seq)
    >>> get_site_complementarity(Seq("ATCG"), Seq("TTAG CC"), 0)
     1
    >>> get_site_complementarity(Seq("ATCG"), Seq("T TAGC C"), 1)
     4
     """
    max_comp = 0
    run_comp = 0
    oligo_compliment = oligo.complement()

    # At each base, check for equality between the oligo's compliment and seq
    for i in range(len(oligo)):
        if oligo_compliment[i] == seq[i + offset]:
            run_comp += 1
        else:
            if max_comp < run_comp:
                max_comp = run_comp
            run_comp = 0
    if max_comp < run_comp:
        max_comp = run_comp
    return max_comp
