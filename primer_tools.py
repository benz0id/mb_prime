from typing import List
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


class PrimerSet:
    """A set of complete MBPrimers and their attributes.
    """

    score: int
    forward_primers: List[MBPrimer]
    reverse_primers: List[MBPrimer]

    def __init__(self, forward_primers: List[MBPrimer],
                 reverse_primers: List[MBPrimer]) -> None:
        """Initialises the Primer set with the given list of primers."""
        self.forward_primers = forward_primers
        self.reverse_primers = reverse_primers
        self.score = 0

    def evaluate_heterogeneity_spacers(self) -> None:
        """Evaluates various attributes of the primers"""
        cmplmnt_lst = []

        for primer in self.forward_primers:
            # Calculate max_comp for any forward primer and forward hetero.
            rev_hetero = primer.get_heterogen_seq().__str__()[::-1]
            cmplmnt_lst.append(get_max_complementarity(Seq(rev_hetero),
                                                       self.forward_primers))
            # No need to reverse hetero seq when calculating binding to reverse
            # primer.
            cmplmnt_lst.append(
                get_max_complementarity(primer.get_heterogen_seq(),
                                        self.reverse_primers))

        # Same process for reverse primers
        for primer in self.reverse_primers:
            rev_hetero = primer.get_heterogen_seq().__str__()[::-1]
            cmplmnt_lst.append(get_max_complementarity(Seq(rev_hetero),
                                                       self.reverse_primers))
            cmplmnt_lst.append(
                get_max_complementarity(primer.get_heterogen_seq(),
                                        self.forward_primers))
