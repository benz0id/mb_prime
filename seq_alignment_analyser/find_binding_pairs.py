import collections
import logging
from pathlib import Path
from typing import List, Union

import heapq
import tqdm

from config_handling.formatting import TargetRegionInfo, AdapterPair, \
    PrimerParams
from seq_alignment_analyser.align import MSA
from seq_alignment_analyser.best_primers import HeteroSeqIterator, \
    BindingPairParams
from seq_alignment_analyser.iterator_manager import BindingIteratorManager
from seq_alignment_analyser.scoring import ScoreBindingPair
from seq_alignment_analyser.sequence_management import PrimerPartsManager, \
    BindingPair, as_binding_pair

log = logging.getLogger('root')

def get_prog_iter(ran: range, do_progbar: bool) -> Union[tqdm.std.tqdm, range]:
    """Returns an iterator over the given range. Displays an error bar iif
    <do_progbar>."""
    if do_progbar:
        return tqdm.tqdm(ran)
    else:
        return ran


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

    _cur_target: The current target
    _cur_msa: The MSA on which the current target lies.

    _bp_heap: Max heap that sorts binding pairs based on their unified scores.

    _num_mt_high: The number of binding sequences with melting temperatures too
    _num_mt_low:  low or high to be considered during any given iteration.

    _num_no_gc:   The number of binding pairs with one or more of their binding
        sequences lacking a GC clamp.

    _total_bps: The total number of binding pairs seen this iteration.

    # Misc.
    _show_prog_bar: Whether to show a progress bar.
    """

    _iterator_manager: BindingIteratorManager
    _parts_manager: PrimerPartsManager
    _scorer: ScoreBindingPair

    _show_prog_bar: bool

    _cur_target: str
    _cur_msa: MSA

    _bp_heap: List[BindingPair]

    _num_mt_high: int
    _num_mt_low: int

    _num_no_gc: int

    _total_bps: int

    def __init__(self, target_sites: List[TargetRegionInfo],
                 adapters: List[AdapterPair],
                 primer_params: PrimerParams, alignments_path: Path,
                 targ_mt: float, max_mt_deviance: float,
                 do_prog_bars: bool = True) -> None:
        """Contructs all required attributes and helpers using the given
        values."""

        self._show_prog_bar = do_prog_bars

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

        # Construct parts manager and iterator manager.
        self._parts_manager = PrimerPartsManager(adapters, msa_to_targets)

        self._iterator_manager = BindingIteratorManager(
            msa_to_targets, self._parts_manager.get_binding_pool_alias(),
            *primer_params)

        self._scorer = ScoreBindingPair(self._parts_manager, targ_mt,
                                        max_mt_deviance)

    def get_best_binding_pairs(self) -> List[BindingPair]:
        """Collects binding pairs that amplify the given targets."""
        pass

    def get_best_binding_pair(self, iterator: HeteroSeqIterator) -> BindingPair:
        """Returns the best binding pair from the given iterator."""
        self._cur_target = iterator.target_name
        self._cur_msa = self._iterator_manager.find_msa(iterator.target_name)

    def update_mt_counters(self, f_seq: str, r_seq: str) -> bool:
        """Updates the melting temp counters. Returns whether the given pair of
        sequences is in range of the other melting temps."""
        self._total_bps += 1
        match self._scorer.is_in_mt_range_seqs(f_seq, r_seq):
            case 'h':
                self._num_mt_high += 1
            case 'y':
                return True
            case 'l':
                self._num_mt_low += 1
        return False

    def store_if_high_conservation(self, param: BindingPairParams, n: int) \
        -> None:
        """If the binding pair specified by <param> has a high conservation,
        store it in the heap. Ensures that self._bp_heap does not exceed <n>
        in length."""
        binding_pair = as_binding_pair(param, self._cur_msa,
                                       self._cur_target)
        for_cons, rev_cons = self._scorer.get_weighted_conservation(
            binding_pair)
        avg_cons = (for_cons + rev_cons) / 2

        # Average of forward and reverse conservation scores.
        binding_pair.set_unified_score(avg_cons)

        # If empty heap or better than worst score.
        add_to_heap = len(self._bp_heap) == 0 or \
            self._bp_heap[0].get_unified_score() < avg_cons

        if add_to_heap:
            heapq.heappush(self._bp_heap, binding_pair)

        # Pop if over size limit.
        if len(self._bp_heap) > n:
            heapq.heappop(self._bp_heap)

    def get_n_most_conserved_valid(self, n: int, iterator: HeteroSeqIterator) \
            -> None:
        """Returns the n most highly conserved valid binding pairs.
        Where conservation is a weighted average of the conservation across the
        forward and reverse binding sequences.

        Binding sequences are valid iff they have a GC in the last two bases and
        a melting temp that's not too deviant from the other primers' melting
        temps."""

        log.info('Beginning iteration through all binding sequences for ' +
                 iterator.target_name)

        # Initialise iteration attributes.
        self._bp_heap = []
        self._num_mt_high = 0
        self._num_mt_low = 0
        self._num_no_gc = 0
        self._total_bps = 0

        num_to_iter = iterator.get_num_pos_primers()

        if self._show_prog_bar:
            print('Finding most conserved primers for' +
                  iterator.target_name + '.')
        for _ in get_prog_iter(range(num_to_iter), self._show_prog_bar):
            f_seq, r_seq = iterator.__next__()

            # Is the melting temp is range of the other seqs?
            if not self.update_mt_counters(f_seq, r_seq):
                continue

            # Is there a GC clamp present?
            if not self._scorer.has_gc_clamp(f_seq, r_seq):
                continue

            # Get conservation score and return.
            param = iterator.get_last_param()
            self.store_if_high_conservation(param, n)

        cons_scores = [bp.get_unified_score() for bp in self._bp_heap]
        avg_con_score = sum(cons_scores) / len(cons_scores)
        max_cons_score = max(cons_scores)

        log.info(''.join(
            [
                'Run Info: \n',
                '\tTotal Number of Binding Pairs Evaluated: ',
                str(self._total_bps), '\n',
                '\tNumber With too High a Tm:        ', str(self._num_mt_high), '\n',
                '\tNumber With too Low a Tm:         ', str(self._num_mt_low), '\n',
                '\tNumber With Good Tm no GC Clamp:  ', str(self._num_no_gc),
                '\n',
                '\tAverage Conservation Amoung ', str(n), ' Best: ',
                str(avg_con_score), '\n',
                '\tMax Conservation Amoung ', str(n), ' Best: ',
                str(max_cons_score), '\n',


            ]
            )
        )





































