from typing import List
import primer3
from hetero_spacer_generator.primer_tools import MBPrimer


def construct_bp_primers(fiveps: List[str], heteros: List[str],
                         bindings: List[str]) -> List[MBPrimer]:
    """Constructs MB primers from the given primer components."""
    assert len(fiveps) == len(heteros) == len(bindings)

    primers = []
    for p_ind in range(len(fiveps)):
        primers.append(MBPrimer(adapter_seq=fiveps[p_ind],
                                binding_seq=bindings[p_ind],
                                index_seq='',
                                heterogen_seq=heteros[p_ind]))
    return primers


def get_all_dgs(seqs: List[str]) -> List[float]:
    """Returns the melting temperatures of all possible structures formed by
    any of seqs."""
    dgs = []
    for i, seq1 in enumerate(seqs):
        for seq2 in seqs[i:]:
            dgs.append(primer3.calcHeterodimer(seq1, seq2).dg)

    return dgs


class PrimerPool:
    """Responsible for the construciton and storage of a collection of complete
    metabarcoding primers.

    _f_primers: The forward primers.
    _r_primers: The reverse primers.
        Where _f_primers[i] and _r_primers[i] target the same site.
    """
    _f_primers: List[MBPrimer]
    _r_primers: List[MBPrimer]

    def __init__(self, f_5p: List[str], f_hetero: List[str],
                 f_binding: List[str], r_5p: List[str], r_hetero: List[str],
                 r_binding: List[str]) -> None:
        """Combines the given components and places them into this primer pool.
        """
        self._f_primers = construct_bp_primers(f_5p, f_hetero, f_binding)
        self._r_primers = construct_bp_primers(r_5p, r_hetero, r_binding)

    def get_all_seqs(self) -> List[MBPrimer]:
        """Returns all the sequences in this pool 5'-3'"""
        return self._f_primers[:] + self._r_primers[:]

