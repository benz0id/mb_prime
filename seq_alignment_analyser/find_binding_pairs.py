import collections
import logging
from pathlib import Path
from typing import List, Union
from copy import deepcopy
import config_handling.command_line_tools as cli
import heapq
import tqdm

from config_handling.formatting import TargetRegionInfo, AdapterPair, \
    PrimerParams, to_range, incl_to_range
from seq_alignment_analyser.align import MSA
from seq_alignment_analyser.best_primers import HeteroSeqIterator, \
    BindingPairParams
from seq_alignment_analyser.iterator_manager import BindingIteratorManager, \
    RESTRICTED, OPEN
from seq_alignment_analyser.scoring import ScoreBindingPair
from seq_alignment_analyser.sequence_management import PrimerPartsManager, \
    BindingPair, as_binding_pair

log = logging.getLogger('root')

# What percent of binding sequences must be within melting temp range.
RETRY_PERCENT = 25
# How much to adjust the length of a binding seq by in order to reach some
# target melting temp.
ADJ_LEN_INCR = 1
# Number of highest conservation primers to keep.
NUM_TO_KEEP = 1000
# The % conservation willing to be lost in order to avoid dimer structures.
MAX_CONS_LOSS = 0.2
# The minimum heap size post iteration.
MIN_HEAP_SIZE = 20


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

    _num_f_mt_high: The number of binding sequences with melting temperatures too
    _num_f_mt_low:  low or high to be considered during any given iteration. (F)

    _num_r_mt_high: The number of binding sequences with melting temperatures too
    _num_r_mt_low:  low or high to be considered during any given iteration. (R)

    _num_no_gc:   The number of binding pairs with one or more of their binding
        sequences lacking a GC clamp.

    _total_bps: The total number of binding pairs seen this iteration.

    # Misc.
    _show_prog_bar: Whether to show a progress bar.

    _mode: Whether to scan all possible lengths at the same time.
    _max_lens: Boolean representation of mode.

    _binding_lens: Allowable binding lengths.


    """

    _iterator_manager: BindingIteratorManager
    _parts_manager: PrimerPartsManager
    _scorer: ScoreBindingPair

    _mode: str
    _show_prog_bar: bool

    _cur_target: str
    _cur_msa: MSA

    _binding_lens: range

    _bp_heap: List[BindingPair]

    _num_f_mt_high: int
    _num_f_mt_low: int
    _num_r_mt_high: int
    _num_r_mt_low: int

    _num_no_gc: int

    _total_bps: int

    def __init__(self, target_sites: List[TargetRegionInfo],
                 adapters: List[AdapterPair],
                 primer_params: PrimerParams, alignments_path: Path,
                 targ_mt: float, max_mt_deviance: float,
                 aln_type: str = 'fasta', do_prog_bars: bool = True,
                 mode: str = RESTRICTED) -> None:
        """Contructs all required attributes and helpers using the given
        values."""

        params = repr(primer_params)[13:-1].split(', ')
        params[2:4] = [', '.join(params[2:4])]
        params[3:5] = [', '.join(params[3:5])]
        params.append('target_melting_temperature=' + str(targ_mt))
        params.append('max_melting_temp_deviance=' + str(max_mt_deviance))
        params = '\n\t' + '\n\t'.join(params)

        targets = [repr(target) for target in target_sites]
        targets = '\n\t' + '\n\t'.join(targets)

        log.info(''.join(
            [
                'Beggining run.', '\nParams:', params, '\nTargets:', targets

            ]
        ))

        self._binding_lens = incl_to_range(primer_params.binding_region_len)

        self._show_prog_bar = do_prog_bars
        self._mode = mode
        if mode == RESTRICTED:
            self._max_lens = False
        else:
            self._max_lens = True
        # Construct attributes required for helpers.
        msa_to_targets = {}
        msa_name_to_msa = {}

        for target in target_sites:
            msa_name = target.aln_filename

            # Parse and store MSAs.
            if msa_name not in msa_name_to_msa.keys():
                msa = MSA(alignments_path / msa_name, filetype=aln_type)
                msa_name_to_msa[msa_name] = msa
                msa_to_targets[msa] = []
            # Multiple targets on one MSA
            else:
                msa = msa_name_to_msa[msa_name]

            msa_to_targets[msa].append(target)

        # Construct parts manager and iterator manager.
        self._parts_manager = PrimerPartsManager(adapters, msa_to_targets)

        self._iterator_manager = BindingIteratorManager(
            msa_to_targets, self._parts_manager.get_binding_pool_alias(), mode,
            *primer_params)

        self._scorer = ScoreBindingPair(self._parts_manager, targ_mt,
                                        max_mt_deviance)

    def get_best_binding_pairs(self) -> List[BindingPair]:
        """Collects binding pairs that amplify the given targets."""
        best_bps = []

        for iterator in self._iterator_manager:
            bp = self.get_best_binding_pair(iterator)
            self._scorer.add_bp(bp)
            self._parts_manager.add_bp(bp)
            best_bps.append(bp)
        log.info('Final Binding Pairs:\n\t' +
                 '\n\t'.join([repr(bp) for bp in best_bps]))

        return best_bps

    def get_best_binding_pair(self, iterator: HeteroSeqIterator) -> BindingPair:
        """Returns the best binding pair from the given iterator."""
        self._cur_target = iterator.target_name
        self._cur_msa = self._iterator_manager.find_msa(iterator.target_name)

        self.get_n_most_conserved_valid(NUM_TO_KEEP, iterator)

        # See whether a significant number of binding pairs were missed as a
        # consequence of melting temp

        bp_heap_cache = []
        while self.too_many_missed(iterator) and not self._max_lens:
            if not self._adjust_iterator(iterator):
                break
            self.get_n_most_conserved_valid(NUM_TO_KEEP, iterator)
            bp_heap_cache.append(deepcopy(self._bp_heap))

        # If the last heap is empty, try using the previous heaps.
        while not self._bp_heap and bp_heap_cache:
            self._bp_heap = bp_heap_cache.pop(-1)

        if not self._bp_heap:
            if not self._max_lens:
                log.warning(
                    'Failed to find a valid binding pair. Trying again with all'
                    ' lengths at once. Try using less restrictive parameters.')
                self._iterator_manager.max_lens(iterator)
                self._max_lens = True
                bp = self.get_best_binding_pair(iterator)
                self._max_lens = False
                return bp

            raise ValueError(
                'Failed to find binding regions within melting temp range '
                'and/or with a GC clamp present. Try entering more permissive '
                'values, or adjusting target melting temp. See logs for more '
                'detailed description of failure.')

        self._bp_heap = sorted(self._bp_heap)
        self._bp_heap.reverse()
        for bp in self._bp_heap:
            # Higher is better.
            dimer_score = self._scorer.get_worst_dgs_average(bp)
            bp.set_dimer_score(dimer_score)

        if not self._max_lens:
            log.info('Desired melting temp distribution met.')

        max_cons = self._bp_heap[0].get_conservation_score()
        i = 0
        while self._bp_heap[i].get_conservation_score() + MAX_CONS_LOSS > \
            max_cons:
            i += 1

        options = [repr(bp) for bp in self._bp_heap[:i + 1]]
        options = '\n\t' + '\n\t'.join(options)

        log.info(''.join(
            [
                'Selecting lowest dimerisation potential from most highly '
                'conserved binding pairs. Pairs within ', str(MAX_CONS_LOSS),
                '% of the most conserved pair include:', options
            ]
        ))

        max_dg_bp = max(self._bp_heap[:i + 1], key=BindingPair.get_dimer_score)

        log.info('\nFinal Selection for ' + self._cur_target + ': \n\t' +
                 repr(max_dg_bp))

        return max_dg_bp

    def update_mt_counters(self, f_seq: str, r_seq: str) -> bool:
        """Updates the melting temp counters. Returns whether the given pair of
        sequences is in range of the other melting temps."""
        # Unpack melting temp attributes.
        in_mt_range, f_high, f_low, r_high, r_low = \
            self._scorer.is_in_mt_range_seqs(f_seq, r_seq)

        if in_mt_range:
            return True

        # Update counters and return false in case of error.
        self._num_f_mt_high += f_high
        self._num_f_mt_low += f_low
        self._num_r_mt_high += r_high
        self._num_r_mt_low += r_low

        return False

    def too_many_missed(self, iterator: HeteroSeqIterator) -> bool:
        """Returns whether a significant number of primer pairs were missed due
        to the produced melting temps being too high or low."""
        cur_f_len, cur_r_len = self._iterator_manager. \
            get_iterator_primer_size(iterator)

        brl = self._iterator_manager.binding_region_len

        # Have a significant portion of the primers been ignored due to melting
        # temp? Have we already reached a length threshold?
        missed_f = self._num_f_mt_high - self._num_f_mt_low
        if abs(missed_f) / self._total_bps * 100 > RETRY_PERCENT:
            f_len_lim_reached = cur_f_len == min(brl) or cur_f_len == max(brl)
            if not f_len_lim_reached:
                return True

        missed_r = self._num_r_mt_high - self._num_r_mt_low
        if abs(missed_r) / self._total_bps * 100 > RETRY_PERCENT:
            r_len_lim_reached = cur_r_len == min(brl) or cur_r_len == max(brl)
            if not r_len_lim_reached:
                return True

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

        # If empty heap or better than worst score.
        add_to_heap = len(self._bp_heap) < MIN_HEAP_SIZE or \
            self._bp_heap[0].get_unified_score() < avg_cons

        if add_to_heap:
            binding_pair.set_unified_score(avg_cons)
            binding_pair.set_conservation_score(avg_cons)
            heapq.heappush(self._bp_heap, binding_pair)

        # Pop if over size limit.
        if len(self._bp_heap) > n:
            heapq.heappop(self._bp_heap)

    def _adjust_iterator(self, iterator: HeteroSeqIterator) \
            -> bool:
        """Adjusts the iterator to include binding sequences closer to the
        target melting temp."""
        adj_f = False
        f_high = False

        adj_r = False
        r_high = False

        cur_f_len, cur_r_len = self._iterator_manager.\
            get_iterator_primer_size(iterator)

        missed_f = self._num_f_mt_high - self._num_f_mt_low
        if abs(missed_f) / self._total_bps * 100 > RETRY_PERCENT:
            adj_f = True
            f_high = missed_f > 0

        missed_r = self._num_r_mt_high - self._num_r_mt_low
        if abs(missed_r) / self._total_bps * 100 > RETRY_PERCENT:
            adj_r = True
            r_high = missed_r > 0

        def get_change(do_change: bool, do_decrease: bool) -> int:
            if not do_change:
                return 0

            direction = (1, -1)[do_decrease]
            return direction * ADJ_LEN_INCR

        f_len = cur_f_len + get_change(adj_f, f_high)
        r_len = cur_r_len + get_change(adj_r, r_high)

        # Make sure that adjustments do not exceed maximum lengths.
        if f_len not in self._binding_lens:
            f_len = cur_f_len

        if r_len not in self._binding_lens:
            r_len = cur_r_len

        if f_len == cur_f_len and r_len == cur_r_len:
            return False

        log.info(''.join(
            [
                'Retrying run with modified lengths.\n',
                '\tf_len: ', str(cur_f_len), ' -> ', str(f_len), '\n',
                '\tr_len: ', str(cur_r_len), ' -> ', str(r_len),

            ]
        ))

        self._iterator_manager.update_iterator_primer_size(iterator, f_len,
                                                           r_len)
        return True

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
        self._num_f_mt_high = 0
        self._num_f_mt_low = 0
        self._num_r_mt_high = 0
        self._num_r_mt_low = 0
        self._num_no_gc = 0
        self._total_bps = 0

        num_to_iter = iterator.get_num_pos_primers()
        if num_to_iter == 0:
            log.critical('Unable to find binding sequences given the specified '
                         'parameters. Try entering more permissive values.')
            cli.eprint('Unable to find binding sequences given the specified '
                         'parameters. Try entering more permissive values.')
            quit(1)
        num_good = 0

        for _ in get_prog_iter(range(num_to_iter), self._show_prog_bar):
            f_seq, r_seq = iterator.__next__()
            self._total_bps += 1

            # Is the melting temp is range of the other seqs?
            if not self.update_mt_counters(f_seq, r_seq):
                continue

            # Is there a GC clamp present?
            if not self._scorer.has_gc_clamp(f_seq, r_seq):
                self._num_no_gc += 1
                continue

            num_good += 1

            # Get conservation score and return.
            param = iterator.get_last_param()
            self.store_if_high_conservation(param, n)

        try:
            iterator.__next__()
        except StopIteration:
            pass
        else:
            for _ in iterator:
                self._total_bps += 1
            log.debug('Iterator Failure: \n\t' + repr(iterator) + '\n\t' +
                      'Iterator Reports: ' + str(num_to_iter) +
                      '\n\t Actual number: ' + str(self._total_bps))
            assert False

        cons_scores = [bp.get_unified_score() for bp in self._bp_heap]
        if len(cons_scores) > 0:
            avg_con_score = sum(cons_scores) / len(cons_scores)
            max_cons_score = max(cons_scores)
        else:
            avg_con_score = 0
            max_cons_score = 0

        log.info(''.join(
            [
                'Run Info: \n',
                '\tTotal Number of Binding Pairs Evaluated: ',
                str(self._total_bps), '\n',
                '\tNum Good Pairs Found :  ', str(num_good),
                '\n',
                '\tNum Forward With too High a Tm:  ', str(self._num_f_mt_high),
                '\n',
                '\tNum Forward With too Low a Tm:   ', str(self._num_f_mt_low),
                '\n',
                '\tNum Reverse With too High a Tm:  ', str(self._num_r_mt_high),
                '\n',
                '\tNum Reverse With too Low a Tm:   ', str(self._num_r_mt_low),
                '\n',
                '\tNumber With Good Tm no GC Clamp: ', str(self._num_no_gc),
                '\n',
                '\tAverage Conservation Amoung ', str(n), ' Best: ',
                str(avg_con_score), '\n',
                '\tMax Conservation Amoung ', str(n), ' Best: ',
                str(max_cons_score), '\n',


            ]
            )
        )

            









































