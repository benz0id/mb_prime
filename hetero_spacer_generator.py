from typing import List, Tuple, Dict, Any, Callable
from Bio.Seq import Seq
from primer_tools import MBPrimer, MBPrimerBuilder, PrimerSet, HalfSet, \
    eval_consecutive_complementarity, eval_total_complementarity, \
    evaluate_heterogen_binding_cross, remove_highest_scores
from presenters import Presenter, ConsolePresenter
from abc import ABC, abstractmethod
from sequence_tools import get_max_complementarity
import random
from defaults import NUM_SPACERS, MAX_SPACER_LENGTH, NUM_HETERO, \
    NUMBER_TO_CROSS_COMPARE, INITIAL_PRIMER_SET_SIZE, \
    GET_SMALLEST_TOTAL_LEN_DEFAULT, GET_SMALLEST_OF_ANY_SPACER_DEFAULT
# TODO better way to do this


class IncompatibleSpacerError(Exception):
    """Thrown when an incompatible spacer is used during some process."""


# Begin methods for sorting spacers
def get_smallest_total_len_list(spacers: List[Tuple[int, int, int, int]]) \
        -> List[int]:
    """Returns a list mapping the index of a spacer combo to the combined length
    of all spacers in that combo."""
    spacer_combined_length = []
    for spacer in spacers:
        spacer_combined_length.append(sum(spacer))
    return spacer_combined_length


def get_smallest_of_any_spacer_list(spacers: List[Tuple[int, int, int, int]]) \
        -> List[int]:
    """Returns a list mapping the index of a spacer combo to the maximum length
    of a spacer in that combo."""
    spacer_max_len = []
    for spacer in spacers:
        spacer_max_len.append(max(spacer))
    return spacer_max_len


class SpacerAlignmentGen:
    """_max_spacer_length:
            The maximum length of a spacer produced.
    _num_hetero:
            The length of the heterogeneity that should be ensured across.
    _presenter:
            The class that handles the visualisation of data produced.
    _criteria:
            Functions score each spacer combo and return a list of those scores.
            Lower scores are always better.
    _criterion_to_weight:
            A dictionary mapping a criterion to it's relative weight (how
            important the scores returned by that function are when sorting
            spacer combos)."""

    _max_spacer_length: int
    _num_hetero: int
    _presenter: Presenter
    _criteria: List[Callable[[List[Tuple[int, int, int, int]]], List[int]]]
    _criterion_to_weight: \
        Dict[Callable[[List[Tuple[int, int, int, int]]], List[int]], int]

    def _build_criteria(self) -> None:
        """Sets self._criteria to some list of default criteria and
        self._criterion_to_weight to the default weights for those criteria."""

        # Add new criteria here.
        self._criterion_to_weight = {
            get_smallest_of_any_spacer_list: GET_SMALLEST_OF_ANY_SPACER_DEFAULT,
            get_smallest_total_len_list: GET_SMALLEST_TOTAL_LEN_DEFAULT
        }

        self._criteria = list(self._criterion_to_weight.keys())

    def __init__(self, max_spacer_length: int, num_hetero: int,
                 presenter: Presenter = ConsolePresenter()) -> None:
        """Initialises the attributes to the values specified.
        Precondition:
            max_spacer_length >= num_hetero"""
        self._presenter = presenter
        self._max_spacer_length = max_spacer_length
        self._num_hetero = num_hetero
        self._build_criteria()

    def get_all_spacer_combos(self, seq: Seq) \
            -> List[Tuple[int, int, int, int]]:
        """Given a <seq>, will provide a set of alignments of that sequence (
        produced by shifting it to the right) that ensures nucleotide
        diversity across the first <self.num_hetero> bases. Returns a list of
        tuples, each of which contain unique combinations of valid heterogeneity
        spacer lengths. Spacer tuples are ordered from smallest spacer to
        greatest."""

        valid_spacer_combos = []
        spacers = []
        self._get_all_compatible_spacers(seq, spacers, NUM_SPACERS,
                                         spacer_combo_list=valid_spacer_combos)
        return valid_spacer_combos

    def _get_all_compatible_spacers(self, seq: Seq, spacers: List[int],
                                    target_length: int,
                                    spacer_combo_list: List[
                                        Tuple[int, int, int, int]] = None) \
            -> None:
        """Will find all valid spacer combos for <seq> that include <spacers>.
        Will add all valid spacers with a length equal to <target_length> to
        <spacer_combo_list>."""

        if len(spacers) == target_length:
            if spacer_combo_list is not None:
                spacer_combo_list.append(tuple(spacers))
            return None

        for spacer in self._get_compatible_spacers(seq, spacers):
            # The set of spacers to be tested for heterogeneity across
            # self.num_hetero
            trial_spacers = spacers.copy()
            trial_spacers.append(spacer)
            self._get_all_compatible_spacers(seq, trial_spacers, target_length
                                             ,spacer_combo_list)
        return None

    def _get_compatible_spacers(self, seq: Seq,
                                seqs_spacers: List[int]) -> tuple[int, ...]:
        """Returns a tuple containing all spacers < <self.num_hetero> such
        that for all j + spacer < <self.num_hetero>, seq[j + spacer] != any
        seqs[i][j]. Returns an empty tuple if no such spacers exist."""

        valid_spacers = []
        # Assume that there exist no valid spacers less than the greatest
        # already used
        if len(seqs_spacers) > 0:
            start_point = seqs_spacers[-1]
        else:
            start_point = 0

        # Iterate though all possible spacers - add to valid spacers if it is
        # compatible
        for i in range(start_point, self._max_spacer_length + 1):
            if self._is_compatible_spacer(seq, seqs_spacers, i):
                valid_spacers.append(i)

        return tuple(valid_spacers)

    def _is_compatible_spacer(self, seq: Seq, spacers: List[int],
                              spacer: int, ) -> bool:
        """Returns whether <seq> shifted by <spacer> has any matching bases
        in the same position as any sequence in <seqs> - each shifted by their
        respective spacers.

        Preconditions: Every item in spacers is less than spacer.
        """
        # Checking for uniqueness at all bases
        for i in range(0, self._num_hetero - spacer):

            # Checking each of the seqs for matching bases.
            for j in range(0, len(spacers)):
                try:
                    original_spacer_base = seq[i + spacer - spacers[j]]
                    new_spacer_seq_base = seq[i]
                    if original_spacer_base == new_spacer_seq_base:
                        return False
                    # Sequence is shorter than heterogeneity region, spacer
                    # fails to push the sequence out of the heterogeneity
                    # region.
                except IndexError:
                    return False

        return True

    def sort_spacer_combos(self,
                           spacer_combos: List[Tuple[int, int, int, int]]
                           ) -> None:
        """Returns a list of spacer combos sorted according to <self.criteria>.
        """
        # Instantiate a list of scores with the lowest possible scores. This
        # will store the scores for each spacer in <spacer_combos>
        scores = []
        for i in range(len(spacer_combos)):
            scores.append(1)

        for i in range(len(self._criteria)):
            criterion = self._criteria[i]
            weight = self._criterion_to_weight[criterion]
            criterion_scores = criterion(spacer_combos)
            for j in range(len(spacer_combos)):
                scores[j] = scores[j] * criterion_scores[j] * weight

        co_sort(scores, spacer_combos)

    def set_params(self, max_spacer_length: int, num_hetero: int) -> None:
        """Sets the construction parameters to the given values."""
        self._num_hetero = num_hetero
        self._max_spacer_length = max_spacer_length


class HeteroSpacerGen(ABC):
    """_max_spacer_length:
            The maximum length of a spacer produced.
    _num_hetero:
            The length of the heterogeneity that should be ensured across.
    _presenter:
            The class that handles the visualisation of data produced."""
    _max_spacer_length: int
    _num_hetero: int
    _presenter: Presenter

    def __init__(self, max_spacer_length: int, num_hetero: int = NUM_HETERO,
                 presenter: Presenter = ConsolePresenter()) -> None:
        """Initialises the attributes to the values specified."""
        self._presenter = presenter
        self._max_spacer_length = max_spacer_length
        self._num_hetero = num_hetero

    @abstractmethod
    def get_hetero_seqs(self, incomplete_forward_primer: MBPrimerBuilder,
                        incomplete_reverse_primer: MBPrimerBuilder,
                        forward_spacers: Tuple[int, int, int, int],
                        reverse_spacers: Tuple[int, int, int, int],
                        num_to_return: int) -> List[PrimerSet]:
        """Generates <num_to_return> PrimerSets, each containing four complete
        <incomplete_forward_primer> and <incomplete_reverse_primer>, with
        heterogeneity spacers of some lengths specified in <forward_spacers> and
        <reverse_spacers> respectively."""
        pass

    def set_params(self, max_spacer_length: int, num_hetero: int) -> None:
        """Sets the construction parameters to the given values."""
        self._num_hetero_forward = num_hetero
        self._max_spacer_length_forward = max_spacer_length


def add(dic: Dict[Any, List[Any]], key: Any, item: Any) -> None:
    """Adds <item> to the <dic> with <key> if <key> already maps to a list, else
     creates a list at <key> and adds <item> to it."""
    if key in dic.keys():
        dic[key].append(item)
    else:
        dic[key] = [item]


def co_sort(to_sort: List[int], to_follow: List[Any], reverse: bool = False) \
        -> None:
    """Sorts to_sort least to greatest (greatest to least if reverse is true).
    Whenever an item in to_sort is moved index i to index j, the item in
    to_follow at i is moved to j.

    Precondition:
            len(to_follow) >= len(to_sort)"""
    dir = [-1, 1][reverse]
    for i in range(len(to_sort), 1, -1):
        for j in range(i - 1):
            if to_sort[j] * dir < to_sort[j + 1] * dir:
                to_sort[j], to_sort[j + 1] = to_sort[j + 1], to_sort[j]
                to_follow[j], to_follow[j + 1] = to_follow[j + 1], to_follow[j]


def _select_and_set(potential_bases: List[str],
                    unfilled_bases: List[List[int]],
                    column: int, sequence_array: List[List[str]]) -> None:
    """Selects a random base from <potential_bases> and places it in a
    random base position in <sequence_array>[base position][<column>]
    specified in <unfilled_bases>[<column>]. Removes the selected base from
    potential bases."""
    specific_potential_bases = potential_bases.copy()
    base_index = random.choice(unfilled_bases[column])

    # Avoid repeating the last base in the sequence if possible
    previous_base = sequence_array[base_index][column + 1]
    if len(specific_potential_bases) > 1 and \
            previous_base in specific_potential_bases:
        specific_potential_bases.remove(previous_base)

    selected_base = random.choice(specific_potential_bases)
    sequence_array[base_index][column] = selected_base
    unfilled_bases[column].remove(base_index)
    potential_bases.remove(selected_base)


def _gen_sequence_array(binding_seq: Seq,
                        spacers: Tuple[int, int, int, int]) \
        -> List[List[str]]:
    """Returns a representation of the <spacers> and <binding_seq>.

    >>> RandomSpacerGen._gen_sequence_array(Seq('ATCG'), (0, 1, 2, 3))
        [['A', 'T', 'C', 'G'],
        ['',  'A','T', 'C'],
        ['',  '', 'A', 'T'],
        ['',  '',  '', 'A']]"""
    sequence_array = []  # [row][column]
    for i in range(len(spacers)):
        sequence_array.append([])
        # Fill hetero spacer region with empty slots
        for _ in range(spacers[i]):
            sequence_array[i].append('')  # Indicates an empty slot

        # Occupy remaining region (plus an extra base) with primer
        # Third spacer is always the largest.
        for j in range(0, spacers[3] - spacers[i] + 1):
            if j < len(binding_seq):
                sequence_array[i].append(binding_seq[j])
            else:
                sequence_array[i].append('-')

    return sequence_array


class RandomSpacerGen(HeteroSpacerGen):
    """ Creates randomised sequences to occupy a heterogeneity alignment.
    === Private Attributes ===
    _random_per_align:
            The number of randomly generated heterogeneity sequences to be
            produced for each alignment.
    _number_to_cross_compare:
            The number of forward and reverse spacers to compare when evaluating
            primer hetero-dimer formation.
    _criteria:
            A list of functions or methods that sort a list of potential spacer
             sequences <List[Tuple[Seq]]> for some <MBPrimerBuilder> according
             to some metric, and leave some num <int> of the best scoring sets
             of spacer sequences in the list."""
    _random_per_align: int
    _number_to_cross_compare: int
    _criteria: List[Callable[[List[Tuple[Seq]], MBPrimerBuilder, int], None]]

    def __init__(self, max_spacer_length: int, num_hetero: int,
                 presenter: Presenter = ConsolePresenter(),
                 rigour: int = 1) -> None:
        """Initialises the attributes to the values specified. Sets the sampling
        size attributes (self._random_per_align & elf._number_to_cross_compare)
        in accordance with the rigour.

        Note: increases to <rigour> will exponentially increase the runtime of
        self._get_hetero_seqs, but may increase the quality of the heterogeneity
        sequences produced."""
        super().__init__(max_spacer_length, num_hetero, presenter)

        self._random_per_align = INITIAL_PRIMER_SET_SIZE * rigour
        self._number_to_cross_compare = NUMBER_TO_CROSS_COMPARE * rigour
        self._criteria = [
            self._remove_high_dimer_complementarity,
            self._remove_high_consec_complementarity
        ]

    def get_hetero_seqs(self, incomplete_forward_primer: MBPrimerBuilder,
                        incomplete_reverse_primer: MBPrimerBuilder,
                        forward_spacer: Tuple[int, int, int, int],
                        reverse_spacer: Tuple[int, int, int, int],
                        num_to_return: int) -> List[PrimerSet]:
        """Generates <num_to_return> PrimerSets, trying to minimise
        complementarity between the heterogeneity regions in the primers and
        that of other sequences."""
        # Generate semi-random heterogeneity seqs for the for and rev primers.
        forward_spacer_seqs = self._gen_hetero_set(incomplete_forward_primer,
                                                   forward_spacer)
        reverse_spacer_seqs = self._gen_hetero_set(incomplete_reverse_primer,
                                                   reverse_spacer)

        # Filter out all but prcnt_to_keep% potential heterogeneity seqs.
        prcnt_to_keep = (self._number_to_cross_compare /
                         self._random_per_align * 100)
        self._filter_spacer_sets(incomplete_forward_primer,
                                 incomplete_reverse_primer,
                                 forward_spacer_seqs, reverse_spacer_seqs,
                                 prcnt_to_keep)

        # Find the best combinations of forward and reverse primers with the
        # least spacer region complementarity.
        return self._cross_compare(incomplete_forward_primer,
                                   incomplete_reverse_primer,
                                   forward_spacer_seqs,
                                   reverse_spacer_seqs, num_to_return)

    def _cross_compare(self, incomplete_forward_primer: MBPrimerBuilder,
                       incomplete_reverse_primer: MBPrimerBuilder,
                       forward_spacer_seqs: List[Tuple[Seq]],
                       reverse_spacer_seqs: List[Tuple[Seq]],
                       num_to_return: int) -> List[PrimerSet]:
        """Compares various combinations of spacer sets, raking them from lowest
        inter-set heterogeneity spacer / primer complementarity"""
        forward_sets = []
        reverse_sets = []
        # Convert sets of spacers & incomplete primers into sets of complete
        # primers. Assumes that len(forward_spacer_seqs) == len(
        # reverse_spacer_seqs)
        for i in range(len(forward_spacer_seqs)):
            forward_sets.append(HalfSet(incomplete_forward_primer,
                                        forward_spacer_seqs[i]))
            reverse_sets.append(HalfSet(incomplete_reverse_primer,
                                        reverse_spacer_seqs[i]))

        scores = []
        # Create a matrix containing the scores between any 2 of the forward and
        # reverse primers.
        for f in range(len(forward_sets)):
            scores.append([])
            for r in range(len(reverse_sets)):
                scores[f].append(
                    evaluate_heterogen_binding_cross(forward_sets[f],
                                                     reverse_sets[r]))

        # Best scores sorted from lowest to highest complementarity
        min_scores = []
        # The combination of for and rev primers that produced the scores in
        # min_scores.
        min_combos = []
        # Set the scores to a list of the worst possible scores
        for i in range(num_to_return):
            min_scores.append(999)
            min_combos.append((-1, -1))

        # If a score is better than the worst score in scores, replace that
        # score and combo with the better score and combo. Sort best to worst.
        for f in range(len(scores)):
            for r in range(len(scores[f])):
                if scores[f][r] < min_scores[-1]:
                    min_scores[-1] = scores[f][r]
                    min_combos[-1] = (f, r)
                    co_sort(min_scores, min_combos)
        # min_scores and min_combos now contain the best scores and the
        # combinations of primers that produced them.
        best_sets = []
        for i in range(num_to_return):
            p_set = PrimerSet(forward_sets[min_combos[i][0]],
                              reverse_sets[min_combos[i][1]])
            best_sets.append(p_set)
        return best_sets

    def _filter_spacer_sets(self, incomplete_forward_primer: MBPrimerBuilder,
                            incomplete_reverse_primer: MBPrimerBuilder,
                            forward_spacer_seqs: List[Tuple[Seq]],
                            reverse_spacer_seqs: List[Tuple[Seq]],
                            percent_to_keep: float) -> None:
        """Removes all but <percent_to_keep>% of elements in
        <forward_spacer_seqs> and <reverse_spacer_seqs>. Removes elements based
        on a set of criteria methods in <self.criteria>.
        Note:
            Rounds percentage down, so may return less spacers than expected."""
        # Calculate the amounts by which to decrement the number of spacers each
        # time
        org_len = len(forward_spacer_seqs)
        percent_to_remove_every_iter = ((100 - percent_to_keep) /
                                        len(self._criteria))
        num_to_remove_each_iter = org_len * percent_to_remove_every_iter / 100
        num_to_keep = org_len
        for criterion in self._criteria:
            num_to_keep -= num_to_remove_each_iter
            criterion(forward_spacer_seqs, incomplete_forward_primer,
                      int(num_to_keep))
            criterion(reverse_spacer_seqs, incomplete_reverse_primer,
                      int(num_to_keep))

    def _remove_high_consec_complementarity(self, spacers: List[Tuple[Seq]],
                                            incomplete_primer: MBPrimerBuilder,
                                            num_to_keep: int) -> None:
        """Removes the spacers from <spacers> that have the greatest number of
        consecutive complementary bases."""
        # Map consecutive complementarity of spacers to their index in
        # <spacers>.
        complementarity_to_index = {}
        for i in range(len(spacers)):
            complement = eval_consecutive_complementarity(incomplete_primer,
                                                          spacers[i])
            add(complementarity_to_index, complement, i)

        # Remove the spacers with highest complementarity.
        remove_highest_scores(spacers, complementarity_to_index, num_to_keep)

    def _remove_high_dimer_complementarity(self, spacers: List[Tuple[Seq]],
                                           incomplete_primer: MBPrimerBuilder,
                                           num_to_keep: int) -> None:
        """Removes the spacers from <spacers> with the highest self
        complementarity with <incomplete_primer>. Retains the <num_to_keep>
        spacers with the lowest self complementarity"""
        # Map complementarity of spacers to their index in <spacers>.
        complementarity_to_index = {}
        for i in range(len(spacers)):
            complement = eval_total_complementarity(incomplete_primer,
                                                    spacers[i])
            add(complementarity_to_index, complement, i)

        # Remove the spacers with highest complementarity.
        remove_highest_scores(spacers, complementarity_to_index, num_to_keep)

    def _gen_hetero_set(self, incomplete_primer: MBPrimerBuilder,
                        spacers: Tuple[int, int, int, int]) \
            -> List[Tuple[Seq]]:
        """Generates <self._random_per_align> semi-random heterogeneity spaces.
        Ensures that primers match the lengths specified in <spacers> and ensure
        heterogeneity with the <incomplete_primer>'s binding_seq."""

        random_spacers = []
        binding_seq = incomplete_primer.get_binding_seq()
        for i in range(self._random_per_align):
            random_spacers.append(
                self._gen_heterogeneity_spacers_rand(binding_seq, spacers))

        return random_spacers

    # Begin methods involved in the generation of randomised heterogeneity
    # sequences.

    def _get_vacant_bases(self, spacers: Tuple[int, int, int, int]) -> List[
        List[int]]:
        """Returns a List containing a List of indices of unfilled bases in
        <sequence_array>.
        >>> seq_arr = [['A', 'T', 'C', 'G'],
        ...             ['',  'A','T', 'C'],
        ...             ['',  '', 'A', 'T'],
        ...             ['',  '',  '', 'A']]
        >>> RandomSpacerGen._get_vacant_bases(seq_arr)
         [[1, 2, 3], [2, 3], [3], []]"""
        unfilled_bases = []
        for i in range(spacers[3] + 1):
            unfilled_bases.append([])

        for i in range(NUM_SPACERS):
            for j in range(spacers[i]):
                unfilled_bases[j].append(i)

        return unfilled_bases

    def _get_potential_bases(self, sequence_array: List[List[str]],
                             column: int) -> List[str]:
        """Returns  a list of all the bases that, if inserted into column, would
        maintain heterogeneity along the aligned nucleotide column.
        >>> seq_arr = [['A', 'T', 'C', 'G'],
        ...             ['',  'A','T', 'C'],
        ...             ['',  '', 'A', 'T'],
        ...             ['',  '',  '', 'A']]
        >>> RandomSpacerGen._get_potential_bases(seq_arr, 0)
         ['T', 'C', 'G']
         >>> RandomSpacerGen._get_potential_bases(seq_arr, 2)
         ['G']"""
        possible_bases = ['A', 'T', 'C', 'G']
        for i in range(NUM_SPACERS):
            # The user entered a set of primers that do not guarantee
            # heterogeneity across the heterogeneity region.
            if sequence_array[i][column].isalnum() \
                    and sequence_array[i][column] not in possible_bases:
                raise IncompatibleSpacerError()
            elif sequence_array[i][column].isalnum():
                possible_bases.remove(sequence_array[i][column])
        return possible_bases

    def _gen_heterogeneity_spacers_rand(self, binding_seq: Seq,
                                        spacers: Tuple[int, int, int, int]) -> \
            Tuple[Seq, ...]:
        """Returns a list of randomly generated valid heterogeneity spacers.
        When possible, avoids base repeats on same strand."""

        sequence_array = _gen_sequence_array(binding_seq, spacers)
        unfilled_bases = self._get_vacant_bases(spacers)

        # Iterate through each column in reverse.
        for column in range(len(unfilled_bases) - 2, -1, -1):
            potential_bases = self._get_potential_bases(sequence_array, column)
            # Randomly fill one of the unspecified bases
            while unfilled_bases[column]:
                _select_and_set(potential_bases, unfilled_bases, column,
                                sequence_array)

        # Extract newly generated heterogeneity spacers from the sequence array.
        final_spacers = []
        for i in range(NUM_SPACERS):
            spacer = Seq(str.join('', sequence_array[i][0:spacers[i]]))
            final_spacers.append(spacer)

        return tuple(final_spacers)

    def set_rigour(self, rigour: int) -> None:
        """Increases the sample sizes of the random generation process according to
        <rigour>."""
        if rigour > 1:
            self._random_per_align = INITIAL_PRIMER_SET_SIZE * rigour
            self._number_to_cross_compare = NUMBER_TO_CROSS_COMPARE * rigour
        elif rigour < 1:
            self._random_per_align = int(INITIAL_PRIMER_SET_SIZE /
                                         -rigour)
            self._number_to_cross_compare = int(NUMBER_TO_CROSS_COMPARE /
                                                -rigour)
        elif rigour == 0:
            raise ValueError("Rigour cannot be equal to 0")


def visualise_complete_primers(primer_sets: List[PrimerSet]) \
        -> Dict[int, PrimerSet]:
    """Prints all of <primers> to the presenter. Returns a dict mapping the
    value corresponding to each primer printed to that primer."""
    to_print = ''
    p_set_dict = {}
    for i in range(len(primer_sets)):
        p_set_dict[i + 1] = primer_sets[i]
        to_print += ''.join(["Primer Set #", str(i + 1), '\n'])
        to_print += str(primer_sets[i])
    return p_set_dict


class HeteroGen:
    """Manages the the creation of heterogeneity spacers. Facades
    SpacerAlignmentGen and RandomPrimerGen.

    === Private Attributes ===
    _max_spacer_length:
            The maximum length of a spacer produced.
    _num_hetero:
            The length of the heterogeneity that should be ensured across.
     _presenter:
            The class that handles the visualisation of data produced.
    _alignment_gen:
            Produces alignments of the binding sequences so that heterogeneity
            can be ensured between them.
    _primer_gen:
            Generates primer sequences.
            """
    _max_spacer_length: int
    _num_hetero: int
    _presenter: Presenter
    _alignment_gen: SpacerAlignmentGen
    _primer_gen: RandomSpacerGen

    def __init__(self, max_spacer_length: int = MAX_SPACER_LENGTH,
                 num_hetero: int = NUM_HETERO,
                 presenter: Presenter = ConsolePresenter()) -> None:
        """Initialises the attributes to the values specified."""
        self._presenter = presenter

        self._max_spacer_length = max_spacer_length
        self._num_hetero = num_hetero
        self._alignment_gen = SpacerAlignmentGen(max_spacer_length, num_hetero,
                                                 presenter)
        self._primer_gen = RandomSpacerGen(max_spacer_length, num_hetero,
                                           presenter)

    def set_params(self, max_spacer_length: int, num_hetero: int) -> None:
        """Sets the construction parameters to the given values."""
        self._num_hetero = num_hetero
        self._max_spacer_length = max_spacer_length
        self._primer_gen.set_params(max_spacer_length, num_hetero)
        self._alignment_gen.set_params(max_spacer_length, num_hetero)

    def set_rigour(self, rigour: int) -> None:
        """Sets the rigour to the specified value."""
        self._primer_gen.set_rigour(rigour)

    def get_all_spacer_combos(self, seq: Seq
                              ) -> List[Tuple[int, int, int, int]]:
        """Given a <seq>, will provide a set of alignments of that sequence (
        produced by shifting it to the right) that ensures nucleotide
        diversity across the first <self.num_hetero> bases. Returns a list of
        tuples, each of which contain unique combinations of valid heterogeneity
        spacer lengths. Will sort the alignments in order to provide shorter
        spacers."""
        aligns = self._alignment_gen.get_all_spacer_combos(seq)
        self._alignment_gen.sort_spacer_combos(aligns)
        return aligns

    def get_hetero_seqs(self, incomplete_forward_primer: MBPrimerBuilder,
                        incomplete_reverse_primer: MBPrimerBuilder,
                        forward_spacer: Tuple[int, int, int, int],
                        reverse_spacer: Tuple[int, int, int, int],
                        num_to_return: int) -> List[PrimerSet]:
        """Generates <num_to_return> PrimerSets, trying to minimise
        complementarity between the heterogeneity regions in the primers and
        that of other sequences."""
        return self._primer_gen.get_hetero_seqs(incomplete_forward_primer,
                                                incomplete_reverse_primer,
                                                forward_spacer, reverse_spacer,
                                                num_to_return)

    def visualise_spacer_alignments(self,
                                    spacers: List[Tuple[int, int, int, int]],
                                    seq: Seq) \
            -> Dict[int, Tuple[int, int, int, int]]:
        """Displays a visual representation of all of the <spacers> for the
        given <seq>. Assumes a set of 4 spacers per tuple."""
        spacer_dict = {}
        to_print = ''
        for i in range(0, len(spacers)):
            spacer_tup = spacers[i]
            spacer_dict[i + 1] = spacers[i]
            to_print += 'Spacer #' + str(i + 1) + ' ' + str(spacer_tup) + '\n'

            for j in range(1, min(self._num_hetero, 11)):
                to_print += str(j) + ' '
            to_print += '\n'

            for j in range(0, NUM_SPACERS):
                spacer = spacer_tup[j]

                to_print += '+ ' * spacer
                if spacer != self._max_spacer_length:
                    to_print += ' '.join(
                        str(seq[0:self._num_hetero - spacer])) + ' '

                to_print += ''.join([''.join(
                    str(seq[self._num_hetero - spacer:])),
                    '\n'])

            to_print += '\n'

        self._presenter.print(to_print)
        return spacer_dict

    def visualise_seq_arr(self, seq_arrs: List[List[List[str]]]) -> None:
        """Displays a visual representation of a given <seq_arrs>.

        Precondition: For all valid i, j, and k,
        len(<seq_arr>[i][j][k]) > self._num_hetero

        Note:
            Replaces all '' in seq_arr with 'N' """
        to_print = ''
        for i in range(len(seq_arrs)):
            for j in range(len(seq_arrs[i])):
                for k in range(len(seq_arrs[i][j])):
                    if seq_arrs[i][j][k] == '':
                        seq_arrs[i][j][k] = 'N'

        for i in range(len(seq_arrs)):
            to_print += ''.join(["Sequence array #", str(i + 1), '\n'])

            for j in range(1, min(self._num_hetero, 11)):
                to_print += str(j) + ' '
            to_print += '\n'

            # Add spaces between the bases in the heterogeneity region
            for j in range(NUM_SPACERS):
                to_print += ' '.join(seq_arrs[i][j][0:self._num_hetero]) + ' '
                to_print += ''.join(seq_arrs[i][j][self._num_hetero:]) + '\n'

            to_print += '\n'

        self._presenter.print(to_print)
