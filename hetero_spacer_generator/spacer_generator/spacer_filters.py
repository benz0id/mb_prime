import abc
from abc import abstractmethod
from typing import Callable, List, Tuple, Dict
from ..primer_types import SimpleCriterion, PairWiseCriterionSingle, \
    SpacerSet

from hetero_spacer_generator.defaults import INITIAL_PRIMER_SET_SIZE, \
    NUM_PAIRINGS_TO_COMP, NUM_SPACERS
from hetero_spacer_generator.primer_tools import HalfSet, HeteroSeqTool, \
    MBPrimerBuilder, MBPrimer, \
    PrimerSet, add, eval_consecutive_complementarity, \
    eval_total_complementarity, evaluate_heterogen_binding_cross, \
    get_cross_iteration_pattern, get_n_lowest, get_n_lowest_matrix, \
    get_these_inds, remove_highest_scores

# Intervals for clearing worst performing HalfSets in PairwiseSpacerSorter
FIRST_INTERVAL = 3
REMOVE_HALFSET_INTERVAL = 2

class SpacerSorter(HeteroSeqTool, abc.ABC):
    """A class that sorts spacers based on some parameters.
    _num_pairings_to_compare:
            The number of pairings of spacers to test against each other when
            evaluating pairings."""
    _num_pairings_to_compare: int

    def __init__(self, max_spacer_length: int, num_hetero: int,
                 num_pairings_to_compare: int) -> None:
        super().__init__(max_spacer_length, num_hetero)
        self._num_pairings_to_compare = num_pairings_to_compare

    @abstractmethod
    def filter_and_make_primer_sets(self,
                                    incomplete_forward_primer: MBPrimerBuilder,
                                    incomplete_reverse_primer: MBPrimerBuilder,
                                    forward_spacers: List[SpacerSet],
                                    reverse_spacers: List[SpacerSet],
                                    num_to_return: int) -> List[PrimerSet]:
        pass


def get_all_arrangements(num_states: int, num_positions: int):
    """Gets all possibles arrangements of a list of length <num_positions> with
    each value in the list being a unique number ST 0 =< n < num states
    Returns a List containing tuples of possible combinations where
    len(tuple) == num_positions"""

    possible_states = []
    arrangements = []
    for i in range(num_states):
        possible_states.append(i)
    get_arrangements_recursive(possible_states, arrangements, [], 0,
                               num_positions)
    return arrangements


def get_arrangements_recursive(possible_states: List[int],
                               arrangements: List[Tuple[int]],
                               arrangement: List[int],
                               depth: int, target_depth: int) -> None:
    """Adds one of <possible_states> to <arrangement>, will recurse until
    <target_depth> is reached, adding an <arrangement> of length <target_depth>
    to <arrangements>."""
    if depth == target_depth:
        arrangements.append(tuple(arrangement))
        return

    # Recurse for each possible state
    for i in range(len(possible_states)):
        next_arrangement = arrangement + [possible_states[i]]
        next_possible_states = possible_states[:]
        next_possible_states.pop(i)
        next_depth = depth + 1
        get_arrangements_recursive(next_possible_states, arrangements,
                                   next_arrangement, next_depth, target_depth)
    return


class SortForPairwise(SpacerSorter):
    """Given sets of forward and reverse heterogeneity sequences and the
    incomplete forward and reverse primer, tries to isolate pairings of primers
    that are most likely to ensure a successful PCR.

    _for_primer:
            An incomplete forward primer.
    _rev_primer:
            An incomplete reverse primer.

    _for_seqs:
            A list of heterogeneity spacer sequences meant to pair with
            _for_primer.
    _rev_seqs:
            A list of heterogeneity spacer sequences meant to pair with
            _rev_primer.

    _single_criteria:
            A list of callables that evaluate some aspect of a spacer
            set in conjunction with some incomplete primer. Returns a list of
            scores equal to the number of spacer_sets.
    _single_criteria_weights:
            A dictionary mapping some criterion to its weight.(i.e. its relevant
            importance when evaluating primers.

    _for_single_scores:
            The scores of the forward primer sets.
    _rev_single_scores:
            The scores of the reverse primer sets.

    _pair_criteria:
            A list of callables that rate a list of forward and reverse spacers
            according to some criteria, where the first <SpacerSets> contains
            the forward spacer seqs and the first <MBPrimerBuilder> contains the
            incomplete forward primer. Returns a list of the scores.
    _pair_criteria_weights:
            A dictionary mapping some criterion to its weight.(i.e. its relevant
            importance when evaluating primers.

    _pair_scores:
            A matrix mapping the indices of forward primers (i) and indices of
            reverse primers (j) to the score of that pairing
            (_pair_scores[i][j]).

    _primer_arrangements:
            A list of all potential primer arrangements, where for each
            primer_arrangement in _primer_arrangements and for any set of
            forward and reverse primers forward_primers[i] pairs with
            reverse_primers[primer_arrangement[i]]."""

    _rigour: int

    _for_primer: MBPrimerBuilder
    _rev_primer: MBPrimerBuilder

    _for_seqs: List[SpacerSet]
    _rev_seqs: List[SpacerSet]

    _single_criteria: List[PairWiseCriterionSingle]
    _single_criteria_weights: Dict[PairWiseCriterionSingle, int]

    _for_single_scores: List[int]
    _rev_single_scores: List[int]

    _pair_criteria: List[SimpleCriterion]
    _pair_criteria_weights: Dict[SimpleCriterion, int]

    def __init__(self, max_spacer_length: int, num_hetero: int,
                 rigour: int = 1) -> None:
        """May initialise a functionally empty class. Build method needs to be
        called in order to ensure proper method functionality."""
        super().__init__(max_spacer_length, num_hetero,
                         rigour * NUM_PAIRINGS_TO_COMP)
        self._rigour = rigour
        self._for_primer = MBPrimerBuilder()
        self._rev_primer = MBPrimerBuilder()
        self._for_seqs = []
        self._rev_seqs = []
        self._for_single_scores = []
        self._rev_single_scores = []
        self._pair_criteria = []

    @classmethod
    def build_full(cls, max_spacer_length: int, num_hetero: int,
                   for_primer: MBPrimerBuilder, rev_primer: MBPrimerBuilder,
                   for_seqs: List[SpacerSet], rev_seqs: List[SpacerSet],
                   rigour: int = 1) -> None:
        """Initialises the classes attributes to the given values.
        Precondition:
                len(for_seqs) == len(rev_seqs)"""
        super().__init__(max_spacer_length, num_hetero,
                         rigour * NUM_PAIRINGS_TO_COMP)
        cls._rigour = rigour
        cls._primer_arrangements = get_all_arrangements(4, 4)
        cls._for_primer = for_primer
        cls._rev_primer = rev_primer
        cls._for_seqs = for_seqs
        cls._rev_seqs = rev_seqs
        cls._for_single_scores = []
        cls._rev_single_scores = []
        for i in range(len(for_seqs)):
            cls._rev_single_scores.append(1)
            cls._for_single_scores.append(1)

        cls._pair_criteria = []

    # INSERT NEW CRITERIA HERE
    def build_criteria(self) -> None:
        """Constructs the criteria attributes with some Callables"""
        self._single_criteria = [
            # TODO fill me in
        ]
        self._pair_criteria = [
            # TODO fill me in
        ]

    def evaluate_scores_single(self):
        """Evaluates the scores of both the forward and reverse primers for each
         criterion in <self._single_criteria>"""
        for criterion in self._single_criteria:
            self._get_scores_single(self._for_primer, self._for_seqs, criterion,
                                    self._for_single_scores)
            self._get_scores_single(self._rev_primer, self._rev_seqs, criterion,
                                    self._rev_single_scores)

    def _get_scores_single(self, primer: MBPrimerBuilder,
                           spacer_seqs: List[SpacerSet],
                           criterion: PairWiseCriterionSingle,
                           scores: List[int]) -> None:
        """Updates the score of each SpacerSet in <scores> using the <criterion>
        method."""
        weight = self._single_criteria_weights[criterion]
        for i in range(len(spacer_seqs)):
            scores[i] += criterion(spacer_seqs[i], primer) * weight

    def sort_and_trim(self):
        """Removes worst scoring forward and reverse primers and their scores
        from their respective attributes."""
        best_for_inds, best_for_scores = \
            get_n_lowest(self._for_single_scores,
                         self._num_pairings_to_compare)
        best_rev_inds, best_rev_scores = \
            get_n_lowest(self._rev_single_scores,
                         self._num_pairings_to_compare)

        self._for_seqs, self._for_single_scores = \
            get_these_inds(best_for_inds, self._for_seqs), best_for_scores
        self._rev_seqs, self._rev_single_scores = \
            get_these_inds(best_rev_inds, self._rev_seqs), best_rev_scores


    def filter_and_make_primer_sets(self,
                                    incomplete_forward_primer: MBPrimerBuilder,
                                    incomplete_reverse_primer: MBPrimerBuilder,
                                    forward_spacers: List[SpacerSet],
                                    reverse_spacers: List[SpacerSet],
                                    num_to_return: int) -> List[PrimerSet]:
        """Returns the PrimerSets best suited for pairwise PCR for
        metabarcoding."""
        self.build_full(self._max_spacer_length, self._num_hetero,
                        incomplete_forward_primer, incomplete_reverse_primer,
                        forward_spacers, reverse_spacers, self._rigour)
        self.evaluate_scores_single()
        self.sort_and_trim()


    def filter_primer_sets(self) -> None:

        self.construct_half_sets()
        self.construct_full_set_matrix()
        # Symmetric iteration pattern.
        iter_pttrn = get_cross_iteration_pattern(self._num_pairings_to_compare)

        #
        for i in range(self._num_pairings_to_compare):
            if i - FIRST_INTERVAL >= 0 \
                and i - FIRST_INTERVAL % REMOVE_HALFSET_INTERVAL == 0:
                self.remove_worst_performing()
            for j in range(self._num_pairings_to_compare):
                f, r = iter_pttrn[i][j]




    def construct_half_sets(self) -> Tuple[List[HalfSet], List[HalfSet]]:
        """Returns (forward_halfset, reverse_halfset) constructed from the
        forward and reverse primers and the sets of heterogeneity spacer
        sequences."""
        pass

    def construct_full_set_matrix(self) -> None:
        """Constructs a matrix of PairwisePrimerSets from the lists of forward
        and reverse HalfSets."""
        pass

    def update_averages(self) -> None:
        """Updates the averages of the HalfSets."""

    def remove_worst_performing(self) -> None:
        """Removes the worst performing HalfSets from the list of HalfSets to
        use."""


def remove_high_consec_complementarity(spacers: List[SpacerSet],
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


def remove_high_dimer_complementarity(spacers: List[SpacerSet],
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


def cross_compare(incomplete_forward_primer: MBPrimerBuilder,
                  incomplete_reverse_primer: MBPrimerBuilder,
                  forward_spacer_seqs: List[SpacerSet],
                  reverse_spacer_seqs: List[SpacerSet],
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

    best_sets = []
    min_combos = get_n_lowest_matrix(scores, num_to_return)[0]
    for i in range(num_to_return):
        p_set = PrimerSet(forward_sets[min_combos[i][0]],
                          reverse_sets[min_combos[i][1]])
        best_sets.append(p_set)
    return best_sets


class SortForSimultaneous(SpacerSorter):
    """Sorts spacers according to some parameters. Returns sets of primers that
    would bind minimally when all used simultaneously in PCR.

    === Private Attributes ===
    _number_to_cross_compare:
            The number of forward and reverse spacers to compare when evaluating
            primer hetero-dimer formation.
    _criteria:
            A list of functions or methods that sort a list of potential spacer
             sequences <SpacerSets> for some <MBPrimerBuilder> according
             to some metric, and leave some num <int> of the best scoring sets
             of spacer sequences in the list."""

    _criteria: List[Callable[[List[SpacerSet], MBPrimerBuilder, int], None]]

    def __init__(self, max_spacer_length: int, num_hetero: int,
                 rigour: int = 1) -> None:
        """Initialises the attributes to the values specified. Sets the sampling
        size attributes (self._random_per_align & self._number_to_cross_compare)
        in accordance with the rigour.

        Note: increases to <rigour> will exponentially increase the runtime of
        self._get_hetero_seqs, but may increase the quality of the heterogeneity
        sequences produced."""
        super().__init__(max_spacer_length, num_hetero,
                         NUM_PAIRINGS_TO_COMP * rigour)
        self._criteria = [
            remove_high_dimer_complementarity,
            remove_high_consec_complementarity
        ]

    def filter_and_make_primer_sets(self,
                                    incomplete_forward_primer: MBPrimerBuilder,
                                    incomplete_reverse_primer: MBPrimerBuilder,
                                    forward_spacers: List[SpacerSet],
                                    reverse_spacers: List[SpacerSet],
                                    num_to_return: int) -> List[PrimerSet]:
        """Generates <num_to_return> PrimerSets, trying to minimise
        complementarity between the heterogeneity regions in the primers and
        that of other sequences.
        Note:
            Will mutate forward and reverse spacers."""

        # Filter out all but prcnt_to_keep% potential heterogeneity seqs.
        prcnt_to_keep = (self._num_pairings_to_compare /
                         len(forward_spacers) * 100)
        self._filter_spacer_sets(incomplete_forward_primer,
                                 incomplete_reverse_primer,
                                 forward_spacers, reverse_spacers,
                                 prcnt_to_keep)

        # Find the best combinations of forward and reverse primers with the
        # least spacer region complementarity.
        return cross_compare(incomplete_forward_primer,
                             incomplete_reverse_primer,
                             forward_spacers,
                             reverse_spacers, num_to_return)

    def _filter_spacer_sets(self, incomplete_forward_primer: MBPrimerBuilder,
                            incomplete_reverse_primer: MBPrimerBuilder,
                            forward_spacer_seqs: List[SpacerSet],
                            reverse_spacer_seqs: List[SpacerSet],
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
