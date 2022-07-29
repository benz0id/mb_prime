import collections
from pathlib import Path
from typing import List

from config_handling.formatting import TargetRegionInfo, AdapterPair, \
    PrimerParams
from seq_alignment_analyser.align import MSA
from seq_alignment_analyser.iterator_manager import BindingIteratorManager
from seq_alignment_analyser.scoring import ScoreBindingPair
from seq_alignment_analyser.sequence_management import PrimerPartsManager, \
    BindingPair


class FindBindingPairs:
    """Iterates over several possible primer binding sequences for the given
    set of targets to be multiplexed.

    === Private Attributes ===

    # Helper Classes

    _iterator_manager: Creates and stores iterators that will allow for
        iteration over every possible binding pair for each target.

    _parts_manager: Stores components of the primers and growing selection of
        binding sequences. Allows for retrieval of specific sequence types.

    _scorer: Scores binding pairs based on their thermodynamic properties.

    # Binding pair iteration attributes.

    _bp_heap: Max heap that sorts binding pairs based on their unified scores.

    _num_mt_high: The number of binding sequences with melting temperatures too
    _num_mt_low:  low or high to be considered during any given iteration.

    _total_bps: The total number of binding pairs seen this iteration.
    """

    _iterator_manager: BindingIteratorManager
    _parts_manager: PrimerPartsManager
    _scorer: ScoreBindingPair

    _bp_heap: List[BindingPair]

    _num_mt_high: int
    _num_mt_low: int
    _total_bps: int

    def __init__(self, target_sites: List[TargetRegionInfo],
                 adapters: List[AdapterPair],
                 primer_params: PrimerParams, alignments_path: Path) -> None:
        """Contructs all required attributes and helpers using the given
        values."""

        # Construct attributes required for helpers.
        msa_to_targets = {}
        msa_name_to_msa = {}

        for target in target_sites:
            msa_name = target.aln_filename

            # Parse and store MSAs.
            if msa_name not in msa_name_to_msa.keys():
                msa = MSA(alignments_path / msa_name)
                msa_name_to_msa[msa_name] = msa
                msa_to_targets[msa] = []
            # Multiple targets on one MSA
            else:
                msa = msa_name_to_msa[msa_name]

            msa_to_targets[msa].append(target)

        self._iterator_manager = BindingIteratorManager(*primer_params,
                                                        msa_to_targets)










