from Bio.Seq import Seq

class Primer(Seq):
    """Basic methods common among all primers"""

    def __init__(self, seq):
        super().__init__(seq)


class MB_Primer(Primer):
    """Stores the attributes of a primer used in two step PCR for metabarcoding.

    === Private Attributes ===
    _adapter_seq:
            The adapter region of the primer.
    _index_seq:
            The indexing region.
    _heterogen_seq:
            The heterogeneity sequence of the primer.
    _primer_region:
            The loci specific binding region.
    """
    _adapter_seq: Seq
    _index_seq: Seq
    _heterogen_seq: Seq
    _primer_seq: Seq

    def __init__(self, adapter_seq: str, index_seq: str, heterogen_seq: str,
                 primer_region: str) -> None:
        """Constructs a primer sequence with components in the following order:
        <adapter_seq> <index_seq> <heterogen_seq> <primer_region> with each of
        their sequences read left to right."""
        super().__init__(''.join([adapter_seq, index_seq, heterogen_seq,
                                  primer_region]))

        self._adapter_seq = Seq(adapter_seq)
        self._index_seq = Seq(index_seq)
        self._heterogen_seq = Seq(heterogen_seq)
        self._primer_seq = Seq(primer_region)

    def get_adapter_seq(self) -> Seq:
        """Returns self._adapter_seq."""
        return self._adapter_seq

    def get_index_seq(self) -> Seq:
        """Returns self._index_seq."""
        return self._index_seq

    def get_heterogen_seq(self) -> Seq:
        """Returns self._heterogen_seq."""
        return self._heterogen_seq

    def get_primer_seq(self) -> Seq:
        """Returns self._primer_seq."""
        return self._primer_seq




