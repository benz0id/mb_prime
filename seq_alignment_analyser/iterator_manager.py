import collections
from typing import Dict, List, Tuple

from config_handling.formatting import InclRange, TargetRegionInfo
from seq_alignment_analyser.align import MSA
from seq_alignment_analyser.best_primers import HeteroSeqIterator
from seq_alignment_analyser.sequence_management import BindingPair
import logging
import config_handling.command_line_tools as cli
import config_handling.formatting as fmt

log = logging.getLogger('root')

TargetIndecies = collections.namedtuple(
    'TargetIndecies',
    [
        'name', 'target_start', 'target_end', 'min_forward_ind',
        'max_forward_ind', 'min_reverse_ind', 'max_reverse_ind',
        'designated_start', 'designated_end', 'overlap_ignore'
    ])


def overlap(start1: int, stop1: int, start2: int, stop2: int,
            allowable: int = 0) -> bool:
    """Returns whether there in overlap in the given regions."""
    if stop1 < start2 or stop2 < start1:
        return False

    if stop1 < stop2:
        if start1 < start2:
            if stop1 - allowable < start2:
                return False

    if stop2 + allowable < stop1:
        if start2 < start1:
            if stop1 - allowable < start2:
                return False

    if stop1 + allowable < stop2:
        if start1 < start2:
            if stop2 - allowable < start1:
                return False

    return True


def any_overlap(targ1: TargetIndecies, targ2: TargetIndecies) -> str:
    """Returns whether any level of overlap exists between the two targets.
    'l' indicates that targ2 < targ1, 'h' indicates targ1 <= targ2, and
    '' indicates no overlap or ignored overlap."""
    ans = ''
    upper1_lt_lower2 = targ1.max_reverse_ind < targ2.min_forward_ind
    upper2_lt_lower1 = targ2.max_reverse_ind < targ1.min_forward_ind
    lower1_lt_lower2 = targ1.min_forward_ind < targ2.min_forward_ind

    # targ2 lies after without overlap
    if upper1_lt_lower2:
        return ''

    # targ2 lies before without overlap
    if upper2_lt_lower1:
        return ''

    # targ1 is below targ2, or targ2 is inside targ1
    if lower1_lt_lower2:
        ans = 'h'

    # targ1 us above targ2, or targ1 is inside targ2
    else:
        ans = 'l'

    ignore = targ1.overlap_ignore and targ2.overlap_ignore

    if ans and ignore:
        return ''

    return ans


class BindingIteratorManager:
    """Manages the iterators that provide information to generate new primer
    sequences. Initialises the iterators with values that guarantee success in
    finding binding sequences.Updates iterators to ensure that conflicts don't
    occur with previously selected primers throughout binding sequence isolation
    process.

    primer_pool: Contains all the binding pairs previously selected by the
        program. Necessary to ensure that future iterators produce pairs
        that conflict with previous pairs.

    iterator_queue: A queue of HeteroBindingIterators, sorted by number of
    available binding pairs L-H.

    iterators: An unsorted list of all constructed iterators.

    msa_to_target_indices: Each MSA and the regions being targeted by primers on
        it.

    primer_primer_dist: Min distance between any two primers.
    """

    primer_pool: List[BindingPair]

    iterator_queue: collections.deque[HeteroSeqIterator]

    iterators: List[HeteroSeqIterator]

    msa_to_target_indices: Dict[MSA, List[TargetIndecies]]

    target_region_len: range

    primer_primer_dist: int

    t: str

    def __init__(self, msa_to_targets: Dict[MSA, List[TargetRegionInfo]],
                 primer_pool: List[BindingPair],
                 primer_primer_distance: int, primer_target_distance: int,
                 target_region_len: InclRange, binding_region_len: InclRange,
                 ideal_binding_size: int) -> None:
        """Initialises all iterators 'maximally' i.e. They can select from any
        possible pair that would leave every other target with a valid pair.

        primer_primer_distance: the distance between any two primers.
        primer_target_distance: the distance between any primer and a target
            site.
        target_region_len: the allowable range of lengths for a target region.
        binding_region_len: the allowable range of lengths for a binding region.
        targets: all of the targets that lie on each MSA.
        """
        num_targs = 0
        for msa in msa_to_targets.keys():
            for _ in msa_to_targets[msa]:
                num_targs += 1

        self.t = ''
        self.primer_pool = primer_pool
        self.target_region_len = fmt.incl_to_range(target_region_len)
        self.primer_primer_dist = primer_primer_distance
        self.iterators = []
        self.iterator_queue = collections.deque()
        self.msa_to_target_indices = {}

        log.info(''.join(
            [
                self.t, 'Beginning iterator construction for ', str(num_targs),
                ' targets across ', str(len(msa_to_targets.keys())),
                ' alignments.'
            ]))

        for msa in msa_to_targets.keys():
            self._designate_regions(primer_primer_distance,
                                    primer_target_distance,
                                    target_region_len, binding_region_len,
                                    msa, msa_to_targets[msa])

        self._construct_iterators(ideal_binding_size)

    def __next__(self) -> HeteroSeqIterator:
        """Returns the next HeteroSeqIterator."""
        next_iterator = self.iterator_queue.popleft()
        self._update_iterator_bounds(next_iterator)
        return next_iterator

    def _lip(self) -> None:
        """Increases the log's indentation level by one."""
        self.t += '    '

    def _lid(self) -> None:
        """Decreases the log's indentation level by one."""
        self.t = self.t[:-4]

    def _designate_regions(
            self, primer_primer_distance: int, primer_target_distance: int,
            target_region_len: InclRange, binding_region_len: InclRange,
            msa: MSA, targets: List[TargetRegionInfo]) -> None:
        """Configures binding iterators for all the given <targets>, which lie
        on <msa>."""
        self._lip()
        targ_regions = []
        self.msa_to_target_indices[msa] = targ_regions

        log.info(self.t + 'Designating ' + str(
            len(targets)) + ' target regions on ' +
                 msa.filename)

        for target in targets:
            min_site = min(target.sites)
            max_site = max(target.sites)

            # Defining the region that must be included in the target.
            target_min = min_site - primer_target_distance
            target_max = max_site + primer_target_distance

            # Defining the regions to which primer can bind.
            r_binding_min = target_max + 1
            f_binding_max = target_min - 1

            r_binding_max = f_binding_max + target_region_len.stop + \
                            2 * binding_region_len.stop
            f_binding_min = r_binding_min - target_region_len.stop - \
                            2 * binding_region_len.stop

            # The minimal region required for this to be a valid target.
            minimal_start = f_binding_max - binding_region_len.stop + 1
            minimal_end = r_binding_min + binding_region_len.stop - 1

            # Region that cannot be infringed upon by other primers.
            designated_start = minimal_start - primer_primer_distance
            designated_end = minimal_end + primer_primer_distance

            targ_inds = TargetIndecies(target.name, target_min, target_max,
                                       f_binding_min, f_binding_max,
                                       r_binding_min, r_binding_max,
                                       designated_start, designated_end,
                                       False)

            # Check that the given regions are valid.
            if not (minimal_start > 0 and minimal_end < len(msa)):
                err_str = ''.join(
                    ['Invalid specification for target region received.',
                     target.name, 'has minimal bounds: (', str(minimal_start),
                     ', ', str(minimal_end), '), which lie outside of the '
                                             'alignment. '])
                log.error(err_str + '\n' + repr(targ_inds))
                raise ValueError(err_str)

            if r_binding_min > r_binding_max:
                err_str = ''.join(
                    ['Invalid specification for target region received.',
                     target.name, 'does not have space in target for '
                                  'desired_features.\n', repr(targ_inds)])
                log.error(err_str)
                raise ValueError(err_str)

            targ_regions.append(targ_inds)

        # Test to see if there is overlap.
        for i, targ in enumerate(targ_regions[:-1]):
            for other_targ in targ_regions[i + 1:]:
                if overlap(targ.designated_start, targ.designated_end,
                           other_targ.designated_start,
                           other_targ.designated_end,
                           allowable=primer_primer_distance):
                    err_msg = ''.join([
                        'Not enough space between ', targ.name, ' and ',
                        other_targ.name, '\n Offending Targets:\n\t',
                        repr(targ),
                        '\n\t', repr(other_targ)])
                    logging.error(err_msg)
                    cli.eprint(err_msg)
                    if cli.yes_no_prompt('Ignore?'):
                        targ.overlap_ignore = True
                        other_targ.overlap_ignore = True
                        continue
                    else:
                        raise ValueError(err_msg)
        self._lid()

    def find_target(self, target_name: str) -> TargetIndecies:
        """Returns the requested target."""
        self._lip()

        target = ''
        for msa in self.msa_to_target_indices.keys():
            for targ in self.msa_to_target_indices[msa]:
                if targ.name == target_name:
                    target = targ
                    break

        if not target:
            raise ValueError('Target does not exist.')

        self._lid()
        return target

    def get_target(self, msa: MSA, target_name: str) -> TargetIndecies:
        """Returns the requested target."""
        self._lip()

        target = ''
        for targ in self.msa_to_target_indices[msa]:
            if targ.name == target_name:
                target = targ
                break

        if not target:
            raise ValueError('Target does not exist.')

        self._lid()
        return target

    def _update_iterator_bounds(self, iterator: HeteroSeqIterator) -> None:
        """Removes all binding indices from iterator that have been made invalid
        by previous primer choices."""
        self._lip()

        log.info(''.join(
            [
                self.t, 'Searching iterator for target \'',
                iterator.target_name, '\' for conflicts with chosen primers.',

            ]
        ))

        iterator_target = self.find_target(iterator.target_name)
        iterator_msa = self.find_msa(iterator.target_name)

        for primer_set in self.primer_pool:
            primer_set_target = self.find_target(primer_set.target_name)
            primer_set_msa = self.find_msa(primer_set_target.name)

            # These targets don't lie on the same MSA.
            if primer_set_msa != iterator_msa:
                continue

            # Ignore overlap?
            if iterator_target.overlap_ignore \
                    and primer_set_target.overlap_ignore:
                continue

            # Extract region used by this primer.
            designated_lower = min(0, primer_set.f_5p - self.primer_primer_dist)
            designated_upper = primer_set.r_5p + self.primer_primer_dist

            # Extract region used by iterator.
            iterator_lower, iterator_upper = \
                iterator.get_forward_reverse_bound()

            # Is there a conflict?
            upper_bound_conflict = designated_lower < max(iterator_upper)
            lower_bound_conflict = min(iterator_lower) < designated_upper
            if not (upper_bound_conflict or lower_bound_conflict):
                continue

            new_upper = min(max(iterator_upper), designated_lower - 1)
            new_iterator_upper = range(min(iterator_upper), new_upper + 1)

            new_lower = max(min(iterator_lower), designated_upper + 1)
            new_iterator_lower = range(new_lower, max(iterator_lower) + 1)

            log.info(''.join(
                [
                    self.t, 'Found a conflict: ', primer_set_target.name, '(',
                    str(designated_lower), ', ', str(designated_upper), ').',
                    self.t, 'Adjusting bounds: ', str(iterator_lower), ', ',
                    str(iterator_upper), '  --->   ', str(new_iterator_lower),
                    ', ', str(new_iterator_upper)
                ]
            ))

            if len(new_iterator_upper) == 0 or len(new_iterator_lower) == 0:
                log.critical('Available space for primers reduced to 0 for this'
                             ' iterator.')
                raise RuntimeError('Produced iterator that does not have '
                                   'sufficient space.')

        self._lid()

    def get_iterator_primer_size(self, iterator: HeteroSeqIterator) -> \
            Tuple[int, int]:
        """Updates the given iterators primer size bounds."""
        f_lens, r_lens = iterator.get_lengths()
        return f_lens[0], r_lens[0]

    def update_iterator_primer_size(self, iterator: HeteroSeqIterator,
                                    f_len: int, r_len: int) -> None:
        """Updates the given iterators primer size bounds."""
        it_targ = self.find_target(iterator.target_name)
        allowable_amp_sizes = self._get_amplicon_size(it_targ, f_len, r_len)
        iterator.new_lens([f_len], [r_len], allowable_amp_sizes)

    def find_msa(self, target_name: str) -> MSA:
        """Given a <target_name>, returns that target's MSA."""
        for msa in self.msa_to_target_indices.keys():
            target_names = [targ.name
                            for targ in self.msa_to_target_indices[msa]]
            if target_name in target_names:
                return msa
        raise ValueError('Target not on any alignment.')

    def _get_overlapping_targets(
            self, msa: MSA, target_name: str) -> Tuple[str, str]:
        """Returns the names of targets that overlap the target
        (preceding, following). Returns '' if no such target exists or the
         overlap has been designated to be ignored."""
        self._lip()
        lower = ''
        higher = ''

        target = self.get_target(msa, target_name)

        for other_target in self.msa_to_target_indices[msa]:
            # Skip if this is the same target.
            if other_target.name == target.name:
                continue

            # Get overlap presence and whether to ignore it: l, h or ''
            overlap = any_overlap(target, other_target)

            if not overlap:
                continue

            if overlap == 'l':
                # if another confilict exists, keep the most agregious one.
                if lower:
                    lower = max([lower, other_target],
                                key=lambda a: a.designated_end)
                else:
                    lower = other_target
            else:
                if higher:
                    higher = min([higher, other_target],
                                 key=lambda a: a.designated_start)
                else:
                    higher = other_target
        self._lid()

        def get_name(s) -> str:
            if isinstance(s, str):
                return s
            else:
                return s.name

        return get_name(lower), get_name(higher)

    def _get_allowable_binding_regions(self, msa: MSA, target: TargetIndecies,
                                       ideal_binding_size: int) \
            -> Tuple[range, range]:
        """Gets available 5p forward and reverse binding indices for the given
        target, ensuring that all targets are given adequate room."""
        self._lip()
        lower_targ_name, higher_targ_name = \
            self._get_overlapping_targets(msa, target.name)

        reps = '\n' + self.t + 'Current Target: ' + repr(target)
        higher_targ = None
        lower_targ = None
        if higher_targ_name:
            higher_targ = self.get_target(msa, higher_targ_name)
            reps += '\n' + self.t + "5' Overlapping target: " + \
                    repr(higher_targ)
        if lower_targ_name:
            lower_targ = self.get_target(msa, lower_targ_name)
            reps += '\n' + self.t + "3' Overlapping target: " + repr(lower_targ)
        if not (lower_targ_name or higher_targ_name):
            ov = '.'
        else:
            ov = ', which has overlap with other targets. '

        log_str = ''.join(
            [
                self.t, 'Configuring allowable binding regions for ',
                target.name, ov, reps
            ])
        log.debug(log_str)

        # Default bounds in case of no overlap. More distal bounds are
        # non-inclusive.
        rev_higher_5p_bound = target.max_reverse_ind + 1
        rev_lower_5p_bound = target.min_reverse_ind + ideal_binding_size - 1
        for_higher_5p_bound = target.max_forward_ind - ideal_binding_size + 1
        for_lower_5p_bound = target.min_forward_ind - 1

        # Adjust bounds to ensure no conflict with other target's designated
        # zone.
        if higher_targ_name:
            rev_higher_5p_bound = min(higher_targ.designated_start,
                                      rev_higher_5p_bound)
        if lower_targ_name:
            for_lower_5p_bound = max(lower_targ.designated_end,
                                     for_lower_5p_bound)

        forward_range = range(for_lower_5p_bound + 1, for_higher_5p_bound + 1)
        reverse_range = range(rev_lower_5p_bound, rev_higher_5p_bound)

        log.debug(''.join(
            [
                self.t, 'Isolated binding ranges: Forward ',
                repr(forward_range), ', Reverse: ', repr(reverse_range)
            ]))

        self._lid()
        return forward_range, reverse_range

    def _get_amplicon_size(self, target:TargetIndecies,
                           f_len: int,  r_len: int) -> range:
        """Converts target region to amplicon length."""
        adj = f_len + r_len
        min_len = target.min_reverse_ind - target.max_forward_ind - 1
        min_len = max(min(self.target_region_len), min_len)

        return range(min_len + adj, max(self.target_region_len) + adj + 1)

    def _construct_iterators(self, ideal_binding_size: int) -> None:
        """Constructs iterators to iterate over all primer matching the given
        parameters."""
        self._lip()


        # For every target,
        for msa in self.msa_to_target_indices.keys():
            for target in self.msa_to_target_indices[msa]:

                # Collect allowable binding start inds for a binding sequence
                # with <target_region_length>.
                f_5p_inds, r_5p_inds = \
                    self._get_allowable_binding_regions(msa, target,
                                                        ideal_binding_size)

                # Convert target length to amplicon length. In this context,
                # amplicon does not include 5' seqs.
                allowed_amp_lens = \
                    self._get_amplicon_size(target, ideal_binding_size,
                                            ideal_binding_size)

                # Construct iterator.
                iterator = HeteroSeqIterator(
                    msa.get_consensus(), f_5p_inds, r_5p_inds,
                    [ideal_binding_size], [ideal_binding_size],
                    allowed_amp_lens, target_name=target.name)

                self.iterators.append(iterator)

        log.info(self.t + 'Iterator generation complete.')

        # Sort iterators by number of possible primers.
        sorted_iterators = sorted(self.iterators,
                                  key=HeteroSeqIterator.get_num_pos_primers)

        # Add iterators to queue
        for i, iterator in enumerate(sorted_iterators):
            log.debug(self.t + str(i) + '. ' +
                      str(iterator.get_num_pos_primers()) + ' ' +
                      repr(iterator))
            self.iterator_queue.append(iterator)

        self._lid()
