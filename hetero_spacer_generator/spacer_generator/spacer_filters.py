import abc
from abc import abstractmethod
from typing import Callable, List, Dict, Tuple
from Bio.Seq import Seq

from hetero_spacer_generator.defaults import NUM_PAIRINGS_TO_COMP
from hetero_spacer_generator.primer_tools import HalfSet, HeteroSeqTool, \
    MBPrimerBuilder, MaxInt, PairWiseCriterionSingle, PairwisePrimerSet, \
    PrimerSet, \
    SimpleCriterion, add, \
    eval_consecutive_complementarity, \
    eval_total_complementarity, evaluate_heterogen_binding_cross, \
    get_all_arrangements, get_cross_iteration_pattern, get_n_lowest, \
    get_n_lowest_matrix, \
    get_these_inds, remove_highest_scores, EvalMBPrimer, calculate_score

# Intervals for clearing worst performing HalfSets in PairwiseSpacerSorter
from hetero_spacer_generator.sequence_tools import is_degen

FIRST_INTERVAL = 2
REMOVE_HALFSET_INTERVAL = 2

# A set of four alignments, i.e. the amount by which each primer should be
# shifted to the right in a set of four to ensure heterogeneity.
SpacerAlignment = Tuple[int, int, int, int]

# A list of four sets of forward or reverse heterogeneity spacer sequences.
SpacerSet = Tuple[Seq, Seq, Seq, Seq]



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

    def set_num_pairings_to_compare(self, num_pairings: int) -> None:
        """Sets <self._num_pairings_to_compare> to <num_pairings>."""
        self._num_pairings_to_compare = num_pairings


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
    _pair_criteria_weights_list:
            A list of weights where _pair_criteria_weights_list[i] is the weight
            of _pair_criteria[i]

    _pair_scores:
            A matrix mapping the indices of forward primers (i) and indices of
            reverse primers (j) to the score of that pairing
            (_pair_scores[i][j]).

    _primer_arrangements:
            A list of all potential primer arrangements, where for each
            primer_arrangement in _primer_arrangements and for any set of
            forward and reverse primers forward_primers[i] pairs with
            reverse_primers[primer_arrangement[i]].

    _for_halfsets:
            HalfSets constructed from the best forward SpacerSets.
    _rev_halfsets:
            HalfSets constructed from the best reverse SpacerSets.
    _fullsets:
            Full sets of primers constructed from the forward and reverse
            halfsets.
    """

    _degen: bool

    _primer_evaluator: EvalMBPrimer

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
    _pair_criteria_weights_list: List[int]

    _for_halfsets: List[HalfSet]
    _rev_halfsets: List[HalfSet]

    _fullsets: List[List[PairwisePrimerSet]]

    def __init__(self, max_spacer_length: int, num_hetero: int,
                 num_pairings_to_comp: int = NUM_PAIRINGS_TO_COMP, degen: bool = None) -> None:
        """May initialise a functionally empty class. Build method needs to be
        called in order to ensure proper method functionality."""
        super().__init__(max_spacer_length, num_hetero, num_pairings_to_comp )
        self._degen = degen
        self._rigour = num_pairings_to_comp
        self._for_primer = MBPrimerBuilder()
        self._rev_primer = MBPrimerBuilder()
        self._for_seqs = []
        self._rev_seqs = []
        self._for_single_scores = []
        self._rev_single_scores = []
        self._pair_criteria = []

    def _do_degen_check(self, primer: MBPrimerBuilder) -> None:
        """Checks the sequences to be used in this filter for degeneracy,
        adjusts class attributes accordingly."""
        for seq in primer:
            if is_degen(seq):
                self._degen = True
                return
        self._degen = False
        return

    def _build_full(self, max_spacer_length: int, num_hetero: int,
                    for_primer: MBPrimerBuilder, rev_primer: MBPrimerBuilder,
                    for_seqs: List[SpacerSet], rev_seqs: List[SpacerSet]
                    ) -> None:
        """Initialises the classes attributes to the given values.
        Precondition:
                len(for_seqs) == len(rev_seqs)"""
        if self._degen is not None and not self._degen:
            self._do_degen_check(for_primer)
            self._do_degen_check(rev_primer)
        self._primer_evaluator = EvalMBPrimer(max_spacer_length, num_hetero,
                                              self._degen)
        self._primer_arrangements = get_all_arrangements(4, 4)
        self._for_primer = for_primer
        self._rev_primer = rev_primer
        self._for_seqs = for_seqs
        self._rev_seqs = rev_seqs
        self._for_single_scores = []
        self._rev_single_scores = []
        self._for_halfsets = []
        self._rev_halfsets = []
        for i in range(len(for_seqs)):
            self._rev_single_scores.append(1)
            self._for_single_scores.append(1)

        self._pair_criteria = []
        self._build_criteria()

    # INSERT NEW CRITERIA HERE
    def _build_criteria(self) -> None:
        """Constructs the criteria attributes with some Callables"""
        self._single_criteria = [
            self._primer_evaluator.eval_homo_hetero_spacer_binding_consec,
            self._primer_evaluator.eval_homo_hetero_spacer_binding_total
        ]
        self._single_criteria_weights = {
            self._primer_evaluator.eval_homo_hetero_spacer_binding_consec: 2,
            self._primer_evaluator.eval_homo_hetero_spacer_binding_total: 1
        }
        self._pair_criteria = [
            self._primer_evaluator.eval_hetero_hetero_spacer_binding_consec,
            self._primer_evaluator.eval_hetero_hetero_spacer_binding_total
        ]
        self._pair_criteria_weights = {
            self._primer_evaluator.eval_hetero_hetero_spacer_binding_consec: 2,
            self._primer_evaluator.eval_hetero_hetero_spacer_binding_total: 1
        }
        self._pair_criteria_weights_list = []
        for criterion in self._pair_criteria:
            self._pair_criteria_weights_list.append(
                self._pair_criteria_weights[criterion])

    # Begin methods for finding seqs with lowest homo-dimerisation potential
    # from a pre-existing set of primers.

    def _evaluate_scores_single(self):
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
            set_scores = []
            for seq in spacer_seqs[i]:
                primer.set_heterogen_seq(seq)
                comp_primer = primer.get_mbprimer()
                set_scores.append(criterion(comp_primer))
            scores[i] += calculate_score(set_scores) * weight

    def _sort_and_trim(self):
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
                                    num_to_return: int, degen: bool = None) \
            -> List[PrimerSet]:
        """Returns the PrimerSets best suited for pairwise PCR for
        metabarcoding."""
        self._build_full(self._max_spacer_length, self._num_hetero,
                         incomplete_forward_primer, incomplete_reverse_primer,
                         forward_spacers, reverse_spacers)
        self._evaluate_scores_single()
        self._sort_and_trim()
        self._score_primer_sets()
        return self._get_lowest_scoring_sets(num_to_return)

    def _construct_half_sets(self) -> None:
        """Returns (forward_halfset, reverse_halfset) constructed from the
        forward and reverse primers and the sets of heterogeneity spacer
        sequences."""
        for for_primer_seqs in self._for_seqs:
            self._for_halfsets.append(HalfSet(self._for_primer,
                                              for_primer_seqs))
        for rev_primer_seqs in self._rev_seqs:
            self._rev_halfsets.append(HalfSet(self._rev_primer,
                                              rev_primer_seqs))
        return

    def _construct_full_set_matrix(self) -> None:
        """Constructs a matrix of PairwisePrimerSets from the lists of forward
        and reverse HalfSets."""
        self._fullsets = []
        for f in range(self._num_pairings_to_compare):
            self._fullsets.append([])
            for r in range(self._num_pairings_to_compare):
                self._fullsets[f].append(PairwisePrimerSet(
                    self._for_halfsets[f], self._rev_halfsets[r]))
        return

    def _get_row_average(self, row: int) -> float:
        """Gets the average of the scored fullsets in some <row> of
        <self._fullsets>."""
        sum_fullsets = 0
        num = 0
        for i in range(self._num_pairings_to_compare):
            if self._fullsets[row][i].has_been_scored():
                sum_fullsets += self._fullsets[row][i].get_min_pairing_score()
                num += 1
        return sum_fullsets / num

    def _get_col_average(self, col: int) -> float:
        """Gets the average of the scored fullsets in some column (<col>) of
        <self._fullsets>."""
        sum_fullsets = 0
        num = 0
        for i in range(self._num_pairings_to_compare):
            if self._fullsets[i][col].has_been_scored():
                sum_fullsets += self._fullsets[i][col].get_min_pairing_score()
                num += 1
        return sum_fullsets / num

    def _update_averages(self) -> None:
        """Updates the averages of the HalfSets."""
        for f in range(self._num_pairings_to_compare):
            if self._for_halfsets[f].is_active():
                self._for_halfsets[f].set_avg(self._get_row_average(f))

        for r in range(self._num_pairings_to_compare):
            if self._for_halfsets[r].is_active():
                self._for_halfsets[r].set_avg(self._get_col_average(r))
        return

    def _remove_worst_performing(self) -> None:
        """Removes the worst performing HalfSet from both the forward and
        reverse halfsets. """
        for_scores = []
        rev_scores = []
        for i in range(self._num_pairings_to_compare):
            for_scores.append(self._for_halfsets[i].get_avg())
            rev_scores.append(self._rev_halfsets[i].get_avg())
        worst_for_index = get_n_lowest(for_scores, 1, highest=True)[0][0]
        worst_rev_index = get_n_lowest(rev_scores, 1, highest=True)[0][0]
        self._for_halfsets[worst_for_index].deactivate()
        self._rev_halfsets[worst_rev_index].deactivate()
        return

    def _score_primer_sets(self) -> None:
        """Creates and scores all combinations of
        <self._num_pairings_to_compare> the forward and reverse primer sets."""
        self._construct_half_sets()
        self._construct_full_set_matrix()
        # Symmetric iteration pattern.
        iter_pttrn = get_cross_iteration_pattern(self._num_pairings_to_compare)

        # Iterate through the full set. Is costly to evaluate set with criteria,
        # so continually remove the worst performing ones.
        for i in range(self._num_pairings_to_compare):
            # Remove worst scoring set after a given interval
            if i - FIRST_INTERVAL >= 0 \
                    and (i - FIRST_INTERVAL) % REMOVE_HALFSET_INTERVAL == 0:
                self._update_averages()
                self._remove_worst_performing()
            for j in range(self._num_pairings_to_compare):
                f, r = iter_pttrn[i][j]
                # If sets are still active, do evaluation.
                if self._for_halfsets[f].is_active() \
                        and self._rev_halfsets[r].is_active():
                    self._fullsets[f][r].apply_criteria(
                        self._pair_criteria,
                        self._pair_criteria_weights_list)

    def _get_lowest_scoring_sets(self, num_sets: int) \
            -> List[PairwisePrimerSet]:
        """Constructs a matrix containing the scores of the now evaluated
        PairwisePrimerSets."""
        scores = []
        for f in range(self._num_pairings_to_compare):
            scores.append([])
            for r in range(self._num_pairings_to_compare):
                if not self._for_halfsets[f].is_active() \
                        or not self._rev_halfsets[r].is_active():
                    scores[f].append(MaxInt())
                else:
                    scores[f].append(
                        self._fullsets[f][r].get_min_pairing_score())
        lowest_scoring_inds = get_n_lowest_matrix(scores, num_sets)[0]
        lowest_scoring = []
        for ind in lowest_scoring_inds:
            lowest_scoring.append(self._fullsets[ind[0]][ind[1]])
        return lowest_scoring


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
