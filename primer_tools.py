from typing import Any, List, Tuple, Dict
from sequence_tools import get_max_complementarity
from Bio.Seq import Seq

FORWARD = "forward"
REVERSE = "reverse"


class Primer(Seq):
    """Basic methods common among all primers"""

    def __init__(self, seq: str):
        super().__init__(seq)


class MBPrimer(Primer):
    """Stores the attributes of a primer used in two step PCR for metabarcoding.

    === Private Attributes ===
    _adapter_seq:
            The adapter region of the primer.
    _index_seq:
            The indexing region.
    _heterogen_seq:
            The heterogeneity sequence of the primer.
    _binding_seq:
            The specific binding region.
    _direction:
            Primers are either "forward" or "reverse".

    On Directionality:
            The sequence of the primer will be given in the 5'-3' direction.
    """
    _adapter_seq: Seq
    _index_seq: Seq
    _heterogen_seq: Seq
    _binding_seq: Seq
    _direction: str

    def __init__(self, adapter_seq: str, index_seq: str, heterogen_seq: str,
                 primer_region: str, direction: str = FORWARD) -> None:
        """Constructs a primer sequence with components in the following order:
        <adapter_seq> <index_seq> <heterogen_seq> <primer_region> with each of
        their sequences read left to right."""
        super().__init__(''.join([adapter_seq, index_seq, heterogen_seq,
                                  primer_region]))

        self._adapter_seq = Seq(adapter_seq)
        self._index_seq = Seq(index_seq)
        self._heterogen_seq = Seq(heterogen_seq)
        self._binding_seq = Seq(primer_region)
        self._direction = direction

    def get_adapter_seq(self) -> Seq:
        """Returns adapter_seq."""
        return self._adapter_seq

    def get_index_seq(self) -> Seq:
        """Returns index_seq."""
        return self._index_seq

    def get_heterogen_seq(self) -> Seq:
        """Returns heterogen_seq."""
        return self._heterogen_seq

    def get_binding_seq(self) -> Seq:
        """Returns binding_seq."""
        return self._binding_seq


class MBPrimerBuilder:
    """Stores the attributes of a primer used in two step PCR for metabarcoding.

    === Private Attributes ===
    _adapter_seq:
            The adapter region of the primer.
    _index_seq:
            The indexing region.
    _heterogen_seq:
            The heterogeneity sequence of the primer.
    _binding_seq:
            The locus specific binding region.
    _direction:
            Primers are either "forward" or "reverse".

    On Directionality:
            The sequence of the primer will be given in the 5'-3' direction.
    """
    _adapter_seq: Seq
    _index_seq: Seq
    _heterogen_seq: Seq
    _binding_seq: Seq
    _direction: str

    def set_adapter_seq(self, seq: Seq) -> None:
        """Sets adapter_seq to <seq>."""
        self._adapter_seq = Seq(str(seq))

    def set_index_seq(self, seq: Seq) -> None:
        """Sets index_seq to <seq>."""
        self._index_seq = Seq(str(seq))

    def set_heterogen_seq(self, seq: Seq) -> None:
        """Sets heterogen_seq to <seq>."""
        self._heterogen_seq = Seq(str(seq))

    def set_binding_seq(self, seq: Seq) -> None:
        """Sets binding_seq to <seq>."""
        self._binding_seq = Seq(str(seq))

    def set_direction(self, direction: str) -> None:
        """Sets the direction of this primer to <direction>"""
        self._direction = direction

    def get_adapter_seq(self) -> Seq:
        """Returns adapter_seq."""
        return self._adapter_seq

    def get_index_seq(self) -> Seq:
        """Returns index_seq."""
        return self._index_seq

    def get_heterogen_seq(self) -> Seq:
        """Returns heterogen_seq."""
        return self._heterogen_seq

    def get_binding_seq(self) -> Seq:
        """Returns binding_seq."""
        return self._binding_seq

    def get_direction(self) -> str:
        """Sets the direction of this primer to <direction>"""
        return self._direction

    def get_MBPrimer(self) -> MBPrimer:
        """Returns a completed version of this primer."""
        return MBPrimer(self._adapter_seq.__str__(),
                        self._index_seq.__str__(),
                        self._heterogen_seq.__str__(),
                        self._binding_seq.__str__(),
                        self._direction)


def eval_self_binding(incomplete_primer: MBPrimerBuilder,
                      spacers: Tuple[Seq]) -> int:
    """Returns the max complementarity between any of <spacers> and any
    <incomplete_primer> completed with a spacer in <spacers>."""
    cmplmnt_lst = []
    primers = []
    for spacer in spacers:
        incomplete_primer.set_heterogen_seq(spacer)
        primers.append(incomplete_primer.get_MBPrimer())
    for spacer in spacers:
        # Calculate max_comp for any forward primer and forward hetero.
        rev_hetero = spacer.__str__()[::-1]
        cmplmnt_lst.append(get_max_complementarity(Seq(rev_hetero),
                                                   primers))
        # No need to reverse hetero seq when calculating binding to reverse
        # primer.
    return max(cmplmnt_lst)


def remove_highest_scores(lst: List[Any], scores_dict: Dict[int:List[int]], num_to_keep: int, lowest: bool = False) -> None:
    """Given <scores_dict> which maps scores to a list of indices, removes items
    at all indices in <lst> except for the <num_to_keep> highest scoring items.
    Iff <lowest> is true, then this behaviour is reversed, and the lowest
    scoring items will be kept."""
    best_to_worst = []
    sorted_keys = list(scores_dict.keys())
    sorted_keys.sort()
    # Insert indices such that best_to_worse contains indices in the
    # lowest scoring -> highest scoring direction
    for key in sorted_keys:
        for index in scores_dict[key]:
            best_to_worst.append(index)

    # Score preference is reversed if lower scores are preferred.
    if lowest:
        best_to_worst.reverse()

    # Remove unwanted items at higher indices first.
    ind_to_remove = best_to_worst[num_to_keep: len(best_to_worst)]
    ind_to_remove.sort(reverse=True)

    for index in ind_to_remove:
        lst.pop(index)


class HalfSet:
    """A set of 4 forward or reverse primers.

    === Public Attributes ===
    primers:
            A set of valid Metabarcoding primers"""

    primers: Tuple[MBPrimer]

    def __init__(self, incomplete_primer: MBPrimerBuilder,
                 spacers: Tuple[Seq]) -> None:
        """Constructs a set of primers complete using <incomplete_primer> each
        with one of <spacers>"""
        primers = []
        for spacer in spacers:
            incomplete_primer.set_heterogen_seq(spacer)
            primers.append(incomplete_primer.get_MBPrimer())
        incomplete_primer.set_heterogen_seq(Seq(''))
        self.primers = tuple(primers)

    def __iter__(self):
        """Returns an iterator for the set of primers contained in this set."""
        return iter(self.primers)


def evaluate_heterogen_binding_cross(forward_primers: HalfSet,
                                     reverse_primers: HalfSet) -> int:
    """Returns the max site complementarity between the heterogeneity
    sequences in <forward_primers> and the <reverse_primers> and visa versa. """
    max_comp = 0

    # Calculate the maximum complementarity between any of the reverse and
    # forward primers
    for primer in forward_primers:
        comp = get_max_complementarity(primer.get_heterogen_seq(),
                                       reverse_primers)
        if comp > max_comp:
            max_comp = comp

    # Same process for reverse primers
    for primer in reverse_primers:
        comp = get_max_complementarity(primer.get_heterogen_seq(),
                                       forward_primers)
        if comp > max_comp:
            max_comp = comp

    return max_comp


class PrimerSet:
    """A set of MBPrimers and their attributes.

    === Private Attributes ===

    """

    forward_primers: List[MBPrimer]
    reverse_primers: List[MBPrimer]

    def __init__(self, forward_primers: List[MBPrimer],
                 reverse_primers: List[MBPrimer]) -> None:
        """Initialises the Primer set with the given list of primers."""
        self.forward_primers = forward_primers
        self.reverse_primers = reverse_primers
