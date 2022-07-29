from typing import List
from primer3 import calcTm, calcHeterodimer
from .sequence_management import BindingPair, PrimerPartsManager


class ScoreBindingPair:
    """Responsible for generating unified scores that capture the features of
    binding pairs.

    === Private Attributes ===

    _melting_temps: The melting temperatures of all binding sequences included
        so far.

    _av_mt: Average melting temp amoung selected primers.

    _targ_mt: The target melting temp specified by the user.

    _primer_parts: Stores all components necessary to construct primer sets.

    _max_mt_deviance: The maximum allowable melting temperature
    deviance between any two binding sequences.
    """

    def __init__(self, primer_parts: PrimerPartsManager, targ_mt: float,
                 max_mt_deviance: float) -> None:
        """Initialises the class with the given values."""
        self._primer_parts = primer_parts
        self._targ_mt = targ_mt
        self._max_mt_deviance = max_mt_deviance

        self._melting_temps = []

    def add_bp(self, bp: BindingPair) -> None:
        """Adds the given binding pair to this classes list of considered
        binding pairs."""
        self._melting_temps.append(bp.f_mt)
        self._melting_temps.append(bp.r_mt)

    def get_avg_conservation(self, bp: BindingPair) -> float:
        """Gets the average conservation across both of the given binding
        regions. """
        f_start = bp.f_5p
        f_end = bp.f_5p + bp.f_len

        r_start = bp.r_5p - bp.r_5p + 1
        r_end = bp.r_5p + 1

        total_cons = bp.msa.get_total_conservation(f_start, f_end) + \
            bp.msa.get_total_conservation(r_start, r_end)

        return total_cons / (bp.f_len + bp.r_len)

    def set_mts(self, bp: BindingPair) -> None:
        """Gets and stores the mating temps of the given <bp>"""
        bp.f_mt = calcTm(bp.f_seq)
        bp.r_mt = calcTm(bp.r_seq)

    def get_mt_deviance(self, bp: BindingPair) -> float:
        """Gets the maximal deviance from the melting temps among the already
        selected primers."""
        # Get maximum deviance.
        max_dev = 0
        for temp in self._melting_temps:
            delta_mt = max([temp - bp.r_mt, temp - bp.f_mt], key=abs)
            if abs(delta_mt) > abs(max_dev):
                max_dev = delta_mt

        delta_mt = bp.r_mt - bp.f_mt
        if abs(delta_mt) > abs(max_dev):
            max_dev = delta_mt

        return max_dev

    def is_in_mt_range(self, bp: BindingPair) -> str:
        """Returns whether <bp> has a melting temp compatible with the already
        selected binding pairs. Where:
         'y' -> Compatible melting temp.
         'h' -> Melting temp too high.
         'l' -> Melting temp too low. """
        max_dev = self.get_mt_deviance(bp)
        if abs(max_dev) <= self._max_mt_deviance:
            return 'y'
        elif max_dev > 0:
            return 'l'
        else:
            return 'h'

    def get_bp_dgs(self, bp: BindingPair) -> List[float]:
        """Returns a list containing the delta g's of formation with all other
        binding pairs, as well as homodimers formed by <bp>."""
        seqs = self._primer_parts.get_all_binding_seqs()

        dgs = []
        for seq in seqs:
            dgs.append(calcHeterodimer(bp.f_seq, seq).dg)
            dgs.append(calcHeterodimer(bp.r_seq, seq).dg)

        dgs.append(calcHeterodimer(bp.f_seq, bp.f_seq).dg)
        dgs.append(calcHeterodimer(bp.r_seq, bp.r_seq).dg)
        dgs.append(calcHeterodimer(bp.r_seq, bp.f_seq).dg)
        return dgs

    def get_5p_seqs_dgs(self, bp: BindingPair) -> List[float]:
        """Returns a list containing all the free energies of formation with 5p
        sequences, including self structures."""
        seqs = self._primer_parts.get_all_5p()
        dgs = []
        for seq in seqs:
            dgs.append(calcHeterodimer(bp.f_seq, seq).dg)
            dgs.append(calcHeterodimer(bp.r_seq, seq).dg)
        return dgs

    def set_unified_score(self, bp: BindingPair):
        """Sets the unified score, which is a combination of all of <bp>'s
        primer properties."""



