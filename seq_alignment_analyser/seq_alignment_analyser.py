from typing import Callable, List
import math
from Bio.Seq import Seq

# The number of nearest bases to consider when calculating the average binding
# of a Seq.
from hetero_spacer_generator.sequence_tools import compare_bases_degenerate, \
    is_degen, is_degen_base

LOCAL_AVERAGE_BREADTH = 3

class SeqAlignAnalyser:
    """A class designed to analyse aligned sequences, and recommend MBPrimer
    binding sites accordingly.

    === Private Attributes ===
    _for_adapter, _for_index, _rev_adapter, _rev_index:
            Components of the metabarcoding primers to be created. Note that
            _rev_adapter and _rev_index are 3' - 5'.

    _proposed_for_binding, _proposed_rev_binding:
            Proposed binding sequences for binding to the _seqs. Note that
            _proposed_rev_binding is 3' - 5'

    _seqs:
            The set of aligned sequences to be compared. 5' - 3'.
            All seqs must have identical length.
    _seqs_complement:
            The complements of each seq.
    _seqs_degen_array:
            Stores the locations of all degenerate bases within seqs.
    _seqs_hg_array:
            Stores the outputs of base comparisons between every seq in seqs at
            each base.
    _seqs_num_matching:
            The number of bases that match between any two seqs at each base.
    _relative_hg:
            The relative heterogeneity at each base.
    _seqs_hg_local_average:
            The average heterogeneity of LOCAL_AVERAGE_BREADTH bases about any
            given base.

    _forward_binding_start, _forward_binding_end, _amplicon_start,
    _amplicon_end, _rev_binding_start, _rev_binding_end:
            Indices at which selected regions in the seqs begin. Note
            directionality of _rev_binding region.

    _seqs_updated:
            Stores whether the derived attributes of the _seqs are
            representative of the current seqs. i.e. if _seqs has been changed
            since those values were last calculated.
    _primer_regions_updated:
            Same as _seqs_updated, but in regards to the _primer_region
            attributes and derived attributes.
    _primer_seqs_updated:
            Same as above but for the sequences of the primers.
    """

    _for_adapter: Seq
    _for_index: Seq
    _rev_adapter: Seq
    _rev_index: Seq

    _proposed_for_binding: Seq
    _proposed_rev_binding: Seq

    _seqs: List[Seq]
    _seqs_complement: List[Seq]
    # self._seqs_degen_array[seq_ind][base_ind] -> bool: base is degenerate
    _seqs_degen_array: List[List[bool]]
    _seqs_hg_array: List[List[List[bool]]]
    _seqs_num_matching: List[int]
    _seqs_relative_hg: List[float]
    _seqs_hg_local_average: List[float]

    # All end positions are non-inclusive.
    _forward_binding_start: int
    _forward_binding_end: int
    _amplicon_start: int
    _amplicon_end: int
    _rev_binding_start: int
    _rev_binding_end: int

    _seqs_current: bool
    _primer_regions_current: bool
    _primer_seqs_current: bool

    # self._base_comp(seq1_ind, b1_ind, seq2_ind, b2_ind) -> bool: bases match.
    _base_matching_method: Callable[[int, int, int, int], bool]

    def __init__(self) -> None:
        """Instantiates the class with empty values."""
        self._for_adapter = Seq('')
        self._for_index = Seq('')
        self._rev_adapter = Seq('')
        self._rev_index = Seq('')

        self._proposed_for_binding = Seq('')
        self._proposed_rev_binding = Seq('')

        self._seqs = []
        self._seqs_complement = []
        self._seqs_degen_array = []
        self._seqs_hg_array = []
        self._seqs_num_matching = []
        self._seqs_relative_hg = []
        self._seqs_hg_local_average = []

        self._seqs_current = False
        self._primer_regions_current = False
        self._primer_seqs_current = False

        self._base_matching_method = self._base_matches_generic

    def do_seq_analysis(self) -> None:
        """Analyses all properties of self._seqs and updates all attributes to
        match."""
        self._gen_complements()
        self._do_degen_check()
        self._gen_heterogeneity_array()
        self._gen_num_matching()
        self._gen_seqs_relative_hg()


    def _gen_complements(self) -> None:
        """Generates the complements of <_seqs> in <_seqs_complements>"""
        self._seqs_complement = []
        for seq in self._seqs:
            self._seqs_complement.append(seq.complement())

    def _do_degen_check(self) -> None:
        """Checks to see whether the inputted seqs are degenerate. Will use
        different comp method if this is the case. """
        for seq in self._seqs:
            if is_degen(seq):
                self._gen_degen_array()
                self._base_matching_method = self._base_matches_w_dgn_arr
                return
        self._base_matching_method = self._base_matches_generic

    def _gen_degen_array(self) -> None:
        """Generates the <_seqs_degen_array> using <_seqs>."""
        self._seqs_degen_array = []
        for i in range(len(self._seqs)):
            self._seqs_degen_array.append([])
            for base in self._seqs[i]:
                self._seqs_degen_array[i].append(is_degen_base(base))

    def _base_matches_w_dgn_arr(self, seq1_ind: int, b1_ind: int, seq2_ind: int,
                                b2_ind: int) -> bool:
        """If either of the bases is degenerate, compare using degenerate
        comparison, else use non-degen base comparison method."""
        if self._seqs_degen_array[seq1_ind][b1_ind] or \
            self._seqs_degen_array[seq2_ind][b2_ind]:
            return compare_bases_degenerate(self._seqs[seq1_ind][b1_ind],
                                        self._seqs[seq2_ind][b2_ind])
        else:
            return self._seqs[seq1_ind][b1_ind] == \
                   self._seqs[seq2_ind][b2_ind]

    def _base_matches_generic(self, seq1_ind: int, b1_ind: int, seq2_ind: int,
                              b2_ind: int) -> bool:
        """Returns whether the specified bases match."""
        return self._seqs[seq1_ind][b1_ind] == \
               self._seqs[seq2_ind][b2_ind]

    def _get_min_seq_len(self) -> int:
        """Returns the length of the shortest seq in <_seqs>."""
        seq_lens = []
        for seq in self._seqs:
            seq_lens.append(len(seq))
        return min(seq_lens)

    def _bases_match(self, base_ind: int, seq1_ind: int, seq2_ind: int) -> bool:
        """Uses the heterogeneity map to confirm whether <seq1_ind> and
        <seq2_ind> have matching bases at <base_ind>."""
        if seq1_ind < seq2_ind:
            seq1_ind, seq2_ind = seq2_ind, seq1_ind

        seq1_ind -= 1
        return self._seqs_hg_array[base_ind][seq1_ind][seq2_ind]

    def _gen_heterogeneity_array(self) -> None:
        """Generates the heterogeneity map using seqs."""
        self._seqs_hg_array = []
        for base_ind in range(self._get_min_seq_len()):
            self._seqs_hg_array.append([])
            for seq1_ind in range(1, len(self._seqs)):
                self._seqs_hg_array[base_ind].append([])
                for seq2_ind in range(0, seq1_ind):
                    self._seqs_hg_array[base_ind][seq1_ind - 1].append(
                        self._base_matching_method(seq1_ind, base_ind, seq2_ind,
                                                   base_ind))
        return

    def _gen_num_matching(self) -> None:
        """Generates the <_seqs_num_matching> list from the heterogeneity map.
        """
        self._seqs_num_matching = []
        for base_ind in range(len(self._seqs_hg_array)):
            total = 0
            for match_list in self._seqs_hg_array[base_ind]:
                for match in match_list:
                    total += match
            self._seqs_num_matching.append(total)

    def _gen_seqs_relative_hg(self) -> None:
        """Generates the <_relative_hg> using <_seqs_num_matching>"""
        n = len(self._seqs) - 1
        # The maximum number of pairs of matching bases within seqs.
        max_pot_hg = n * (n - 1) / 2
        self._seqs_relative_hg = []
        for num_match in self._seqs_num_matching:
            # This gives a measure of the heterogeneity at the base position.
            # We see a rapid decrease in (num_match / max_pot_hg) non matching
            # bases are introduced, so we use sqrt() to counteract this.
            self._seqs_relative_hg.append(
                math.sqrt(num_match / max_pot_hg)
            )

        return

    def _gen_seqs_hg_local_average(self) -> None:
        """Generates the local average by calculating some average of the
        LOCAL_AVERAGE_BREADTH nearest bases on either side."""
        self._seqs_hg_local_average = []

    def add_seq(self, seq: Seq) -> None:
        """Adds seq to the aligned seqs."""
        seq_len = len(seq)
        for existing_seq in self._seqs:
            if len(existing_seq) != seq_len:
                raise ValueError("All inputted seqs must be the same length.")
            return
        self._seqs.append(seq)
        self._seqs_current = False

    def clear_seqs(self) -> None:
        """Clears the aligned sequences."""
        self._seqs = []
        self._seqs_current = False

    def set_seqs(self, seqs: List[Seq]) -> None:
        """Sets the aligned seqs to the given seqs."""
        self._seqs = seqs[:]
        for i in range(1, len(seqs)):
            for j in range(i + 1, len(seqs)):
                if len(seqs[i]) != len(seqs[j]):
                    raise ValueError("All inputted seqs must be the same "
                                     "length.")
        self._seqs_current = False

    def _set_for_adapter(self, seq: Seq) -> None:
        """Sets the value of the <for_adapter> sequence."""
        self._primer_seqs_current = False
        self._for_adapter = seq

    def set_rev_index(self, seq: Seq) -> None:
        """Sets the value of the <rev_index> sequence."""
        self._primer_seqs_current = False
        self._rev_index = seq

    def set_for_index(self, seq: Seq) -> None:
        """Sets the value of the <for_index> sequence."""
        self._primer_seqs_current = False
        self._for_index = seq

    def set_rev_adapter(self, seq: Seq) -> None:
        """Sets the value of the <rev_adapter> sequence."""
        self._primer_seqs_current = False
        self._rev_adapter = seq

    def get_for_adapter(self) -> Seq:
        """Returns the value of the <for_adapter> sequence."""
        return self._for_adapter

    def get_for_index(self) -> Seq:
        """Returns the value of the <for_index> sequence."""
        return self._for_index

    def get_rev_index(self) -> Seq:
        """Returns the value of the <rev_index> sequence."""
        return self._rev_index

    def get_rev_adapter(self) -> Seq:
        """Returns the value of the <rev_adapter> sequence."""
        return self._rev_adapter
