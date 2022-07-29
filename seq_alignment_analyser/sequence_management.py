from copy import deepcopy
from typing import Dict, List, Tuple

from seq_alignment_analyser.align import MSA
import config_handling.formatting as fmt


def comp(seq: str) -> str:
    """Returns the complement of the given DNA sequence."""
    new_seq = ''
    for c in seq:
        if c.upper() == 'A':
            new_seq += 'T'
        elif c.upper() == 'T':
            new_seq += 'A'
        elif c.upper() == 'C':
            new_seq += 'G'
        elif c.upper() == 'G':
            new_seq += 'C'
        else:
            raise ValueError('Invalid dna sequence received: ' + seq)
    return new_seq


def rev_comp(seq: str) -> str:
    """Returns the reverse complement fo the given DNA sequence. """
    return comp(seq)[::-1]


class BindingPair:
    """A pair of binding sequences to be incorporated into forward and reverse
    metabarcoding primers.

    === Public Attributes ===

    # Identifiers
    msa: The multiple sequence alignments to which this binding pairs are
        targeted.
    target_name: The name of target captured by this binding pair.

    # Primer Params
    f_seq: The sequence of the forward primer.
    r_seq: The sequence of the reverse primer.

    f_mt: The melting temp of the forward binding sequence.
    r_mt: The melting temp of the reverse binding sequence.

    f_5p: The index of the most 5' base bound by the forward primer on the
        alignment.
    r_5p: The index of the most 5' base bound by the reverse primer on the
        alignment.

    f_len: The length of the forward primer.
    r_len: The length of the reverse primer.

    # Scoring
    unified_score: The score of this primer set as calculated by an external
    scorer.
    """
    msa: MSA
    target_name: str

    f_seq: str
    r_seq: str

    f_5p: int
    r_5p: int

    f_len: int
    r_len: int

    def __init__(self, msa: MSA, target_name: str, f_5p: int, r_5p: int,
                 f_len: int, r_len: int) -> None:
        """Initialises this pair using hte given attributes. Extracts primer
        sequences from the given alignment."""
        self.msa = msa
        self.target_name = target_name

        self.f_5p = f_5p
        self.r_5p = r_5p

        self.f_mt = 0
        self.r_mt = 0

        self.f_len = f_len
        self.r_len = r_len

        self.f_seq = msa.get_consensus()[f_5p: f_5p + f_len]
        self.r_seq = rev_comp(msa.get_consensus()[r_5p - r_len + 1: r_5p + 1])


class PrimerPartsManager:
    """Manages the components to be incorporated into the final primers.
    Contains several methods that allow for the evaluation of dimerisation
    properties between potential binding sequences and the adapters.

    === Private Attributes ===

    _target_to_5p_seq: Maps the name of a target to the 5' seqs to be attached
        to the binding sequences at that target.

    _binding_pair_pool: The set of binding pairs selected so far.
    """
    _fivep_seqs: List[fmt.AdapterPair]
    _binding_pair_pool: List[BindingPair]
    _msa_to_target: Dict[MSA, List[fmt.TargetRegionInfo]]

    def __init__(self, target_to_fivep_seqs: List[fmt.AdapterPair],
                 msa_to_targets: Dict[MSA, List[fmt.TargetRegionInfo]]) -> None:
        """Initialises this class using the given parameters."""
        self._binding_pair_pool = []
        self._fivep_seqs = deepcopy(target_to_fivep_seqs)
        self._msa_to_target = deepcopy(msa_to_targets)



    def add_bp(self, binding_pair: BindingPair) -> None:
        """Creates and returns a binding pair using the given attributes, as
        well as adding it to the growing """
        self._binding_pair_pool.append(binding_pair)

    def get_all_5p(self) -> Tuple[str, ...]:
        """Returns all 5' sequences in the 5'-3' direction in a list."""
        seqs = []
        for adapter_pair in self._fivep_seqs:
            seqs.append(adapter_pair.forward)
            seqs.append(adapter_pair.reverse)

        return tuple(seqs)

    def get_all_binding_seqs(self) -> Tuple[str, ...]:
        """Returns all the currently collected binding sequences 5' - 3'."""
        seqs = []
        for pair in self._binding_pair_pool:
            seqs.append(pair.f_seq)
            seqs.append(pair.r_seq)
        return tuple(seqs)
