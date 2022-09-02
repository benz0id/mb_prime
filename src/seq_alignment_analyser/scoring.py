from typing import List, Tuple

import primer3

from src.hetero_spacer_generator.sequence_tools import get_p3_adapter_float
from math import ceil
from statistics import mean
from .align import MSA
from .sequence_management import BindingPair, PrimerPartsManager

GC_ZONE_SIZE = 2
NUM_TO_COLLECT = 1 / 4

primer3_adapter = get_p3_adapter_float()
calc_dimer = lambda a, b: primer3_adapter.calc_heterodimer_score(a, b)
calc_tm = primer3.calcTm


class ScoreBindingPair:
    """Responsible for generating unified scores that capture the features of
    binding pairs.

    === Private Attributes ===

    _melting_temps: The melting temperatures of all binding sequences included
        so far.

    _targ_mt: The target melting temp specified by the user.

    _primer_parts: Stores all components necessary to construct primer sets.

    _max_mt_deviance: The maximum allowable melting temperature
    deviance between any two binding sequences.

    _weight_dis_map: Contains conservation weight distributions for sequences
        of up to length 100. Maps sequence length to weight distribution tuple.
    """

    def __init__(self, primer_parts: PrimerPartsManager, targ_mt: float,
                 max_mt_deviance: float) -> None:
        """Initialises the class with the given values."""
        self._primer_parts = primer_parts
        self._targ_mt = targ_mt
        self._max_mt_deviance = max_mt_deviance

        self._melting_temps = [targ_mt]

        self._weight_dis_map = {}
        for i in range(5, 101):
            self._weight_dis_map[i] = \
                self.get_conservation_weight_distribution_last_5(i)

    def add_bp(self, bp: BindingPair) -> None:
        """Adds the given binding pair to this classes list of considered
        binding pairs."""
        self._melting_temps.append(calc_tm(bp.get_f_seq()))
        self._melting_temps.append(calc_tm(bp.get_r_seq()))


    def get_conservation_weight_distribution_last_5(self, len: int) \
            -> Tuple[float]:
        """Returns the weight array for a sequence of <len>.

        Gives the last five bases increased weight with more 3' bases given more
        weight.

        first <len> - 5 bases given uniform weight.

        Last 5 bases given weight according to:
            https://www.desmos.com/calculator/zg5iklewtq
        """
        assert len >= 5

        last_5_weight = 0.5
        first_bases_weight = 0.5

        if len <= 15:
            adj_fac = (len - 5) / 10
            last_5_weight *= 2 - adj_fac
            first_bases_weight *= adj_fac

        # Determines the relative weight at each of the last
        distribution_function = lambda x: ((x + 3) / 6) ** (x + 3) + 1
        last_5_dist = [distribution_function(x) for x in range(0, 5)]
        adj = last_5_weight / sum(last_5_dist)
        last_5_dist = [adj * weight for weight in last_5_dist]

        n = (len - 5)
        if n:
            first_n_dist = [first_bases_weight / n] * n
        else:
            first_n_dist = []

        assert 0.9999999 <= sum(first_n_dist + last_5_dist) <= 1.000001

        return tuple(first_n_dist + last_5_dist)

    def get_weighted_conservation(self, bp: BindingPair) -> Tuple[float, float]:
        """Returns the weighted conservation. Gives an elevated weighting to the
        last ten bases."""
        f_cons = self.get_weighted_seq_conservation(
            bp.f_5p, bp.f_len, False, bp.msa
        )
        r_cons = self.get_weighted_seq_conservation(
            bp.r_5p, bp.r_len, True, bp.msa
        )
        return f_cons, r_cons

    def get_weighted_seq_conservation(self, five_p_ind: int, length: int,
                                      rev: bool, msa: MSA) -> float:
        """Returns the weighted conservation of the given region.

        The 10 most 3' bases are given maximal, identical weights, and the
        remaining bases are given exponentially decaying weights based on their
        distance from the 5' end.
        """
        weight_arr = self._weight_dis_map[length]
        msa_inds = msa.get_primer_indices_array(five_p_ind, length, rev)
        cons_arr = msa.get_conservation_across(msa_inds)
        score = 0
        for i, cons in enumerate(cons_arr):
            score += cons * weight_arr[i]
        return score

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

    def get_pair_deviance(self, f_seq: str, r_seq: str) \
            -> Tuple[float, float, float]:
        """Gets the maximal deviance from the melting temps among the already
        selected primers."""

        f_dev = self.get_mt_deviance(f_seq)
        r_dev = self.get_mt_deviance(r_seq)

        pair_dev = abs(calc_tm(f_seq) - calc_tm(r_seq))

        return f_dev, r_dev, pair_dev

    def get_mt_deviance(self, seq: str) -> float:
        """Gets the maximal deviance from the melting temps among the already
        selected primers."""
        r_mt = calc_tm(seq)

        # Get maximum deviance.
        max_dev = 0
        for temp in self._melting_temps:
            delta_mt = r_mt - temp
            if abs(delta_mt) > abs(max_dev):
                max_dev = delta_mt

        return max_dev

    def is_in_mt_range_seqs(self, f_seq: str, r_seq: str) -> \
            Tuple[bool, bool, bool, bool, bool]:
        """Returns whether <bp> has a melting temp compatible with the already
        selected binding pairs. Returns a tuple of booleans containing
        information regarding the melting temps of the pair.

        Attributes stored in tuple in the following order:
        in_mt_range, f_high, f_low, r_high, r_low
        in_mt_range     : Are there melting temp conflicts?
        f_high          : Is the forward seq melting temp higher than avg?
        f_low           : Is the forward seq melting temp lower than avg?
        r_high          : Is the reverse seq melting temp higher than avg?
        r_low           : Is the reverse seq melting temp lower than avg?

        Only the •problem causing* attributes will be true.
        """
        f_dev, r_dev, pair_dev = self.get_pair_deviance(f_seq, r_seq)

        in_mt_range = False
        f_high = False
        f_low = False
        r_high = False
        r_low = False

        def get_tuple() -> Tuple[bool, bool, bool, bool, bool]:
            return in_mt_range, f_high, f_low, r_high, r_low

        def all_lt(iterable, threshold) -> bool:
            for val in iterable:
                if abs(val) > threshold:
                    return False
            return True

        if all_lt((f_dev, r_dev, pair_dev), self._max_mt_deviance):
            in_mt_range = True
            return get_tuple()

        # Is the forward sequence overly deviant?
        if abs(f_dev) > self._max_mt_deviance:
            if f_dev > 0:
                f_high = True
            else:
                f_low = True

        # Is the reverse sequence overly deviant?
        if abs(r_dev) > self._max_mt_deviance:
            if r_dev > 0:
                r_high = True
            else:
                r_low = True

        # Are the pair together overly deviant?
        if pair_dev > f_dev and pair_dev > r_dev and \
                pair_dev > self._max_mt_deviance:
            if f_dev > r_dev:
                f_high = True
                r_low = True
            else:
                f_low = True
                r_high = True

        return get_tuple()

    def get_worst_dgs_average(self, bp: BindingPair) -> float:
        """Calculates the average free energy among the worst structures formed
        by this binding pair."""
        bps_dgs = self.get_bp_dgs(bp)
        five_p_dgs = self.get_5p_seqs_dgs(bp)

        all_structures = bps_dgs + five_p_dgs

        num_structs = len(all_structures)
        num_worst = ceil(num_structs * NUM_TO_COLLECT)

        worst = sorted(all_structures)[:num_worst]
        return mean(worst)

    def get_bp_dgs(self, bp: BindingPair) -> List[float]:
        """Returns a list containing the delta g's of formation with all other
        binding pairs, as well as homodimers formed by <bp>."""
        seqs = self._primer_parts.get_all_binding_seqs()

        dgs = []
        for seq in seqs:
            dgs.append(calc_dimer(bp.get_f_seq(), seq))
            dgs.append(calc_dimer(bp.get_r_seq(), seq))

        dgs.append(calc_dimer(bp.get_f_seq(), bp.get_f_seq()))
        dgs.append(calc_dimer(bp.get_r_seq(), bp.get_r_seq()))
        dgs.append(calc_dimer(bp.get_r_seq(), bp.get_f_seq()))
        return dgs

    def get_5p_seqs_dgs(self, bp: BindingPair) -> List[float]:
        """Returns a list containing all the free energies of formation with 5p
        sequences, including self structures."""
        seqs = self._primer_parts.get_all_5p()
        dgs = []
        for seq in seqs:
            dgs.append(calc_dimer(bp.get_f_seq(), seq))
            dgs.append(calc_dimer(bp.get_r_seq(), seq))
        return dgs

    def has_gc_clamp(self, f_seq: str, r_seq: str) -> bool:
        """Returns whether both sequences have a GC clamp."""
        f_gc_zone = f_seq[-GC_ZONE_SIZE:]
        r_gc_zone = r_seq[-GC_ZONE_SIZE:]

        f_gc = 'G' in f_gc_zone or 'C' in f_gc_zone
        r_gc = 'G' in r_gc_zone or 'C' in r_gc_zone

        return f_gc and r_gc
