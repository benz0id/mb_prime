import random
from abc import abstractmethod
from typing import Dict, List

from hetero_spacer_generator.defaults import INITIAL_PRIMER_SET_SIZE, \
    NUM_PAIRINGS_TO_COMP, NUM_HETERO, NUM_SPACERS
from hetero_spacer_generator.get_random_seqs import gen_hetero_set
from hetero_spacer_generator.spacer_generator.spacer_filters import \
    SortForPairwise, SortForSimultaneous, SpacerAlignment, SpacerSet, \
    SpacerSorter
from hetero_spacer_generator.primer_tools import HeteroSeqTool, \
    MBPrimerBuilder, PrimerSet


class HeteroSpacerGen(HeteroSeqTool):
    """ A class with which a set of valid spacers can be generated given some
    set of valid alignments and incomplete primers"""

    def __init__(self, max_spacer_length: int, num_hetero: int = NUM_HETERO) \
            -> None:
        """Initialises the attributes to the values specified."""
        super().__init__(max_spacer_length, num_hetero)

    @abstractmethod
    def get_hetero_seqs(self, incomplete_forward_primer: MBPrimerBuilder,
                        incomplete_reverse_primer: MBPrimerBuilder,
                        forward_spacers: SpacerAlignment,
                        reverse_spacers: SpacerAlignment,
                        num_to_return: int) -> List[PrimerSet]:
        """Generates <num_to_return> PrimerSets, each containing four complete
        <incomplete_forward_primer> and <incomplete_reverse_primer>, with
        heterogeneity spacers of some lengths specified in <forward_spacers> and
        <reverse_spacers> respectively."""
        pass

    @abstractmethod
    def set_pairwise(self, degen: bool = None) -> None:
        """Sets the primer types to be returned to pairwise primers"""
        pass

    def set_params(self, max_spacer_length: int, num_hetero: int) -> None:
        """Sets the construction parameters to the given values."""
        self._num_hetero_forward = num_hetero
        self._max_spacer_length_forward = max_spacer_length


class RandomSpacerGen(HeteroSpacerGen):
    """ Creates randomised sequences to occupy a heterogeneity alignment.
    === Private Attributes ===
    _number_to_cross_compare:
            The number of forward and reverse spacers to compare when evaluating
            primer hetero-dimer formation.
    _spacer_sorter:
            Helper class that sorts spacers according to some parameters.
    """
    _random_per_align: int
    _num_pairings_to_compare: int
    _spacer_sorter: SpacerSorter
    _rigour: int

    def __init__(self, max_spacer_length: int, num_hetero: int,
                 rigour: int = 1) -> None:
        """Initialises helper classes and attributes."""
        super().__init__(max_spacer_length, num_hetero)
        self._spacer_sorter = SortForSimultaneous(max_spacer_length,
                                                  num_hetero)
        self.set_rigour(rigour)

    def get_random_hetero_seqs(self, incomplete_forward_primer: MBPrimerBuilder,
                               incomplete_reverse_primer: MBPrimerBuilder,
                               forward_spacer: SpacerAlignment,
                               reverse_spacer: SpacerAlignment,
                               num_to_return: int) -> List[PrimerSet]:
        """Generates <num_to_return> PrimerSets, without minimising
        complementarity between the heterogeneity regions in the primers and
        that of other sequences."""
        forward_spacer_seqs = gen_hetero_set(incomplete_forward_primer,
                                             forward_spacer,
                                             self._random_per_align)
        reverse_spacer_seqs = gen_hetero_set(incomplete_reverse_primer,
                                             reverse_spacer,
                                             self._random_per_align)

        # Choose some number of random spacer sequences.
        primer_sets = []
        for i in range(num_to_return):
            forward_primers = []
            for seq in forward_spacer_seqs[i]:
                incomplete_forward_primer.set_heterogen_seq(seq)
                f_primer = incomplete_forward_primer.get_mbprimer()
                forward_primers.append(f_primer)
            reverse_primers = []
            for seq in reverse_spacer_seqs[i]:
                incomplete_reverse_primer.set_heterogen_seq(seq)
                r_primer = incomplete_reverse_primer.get_mbprimer()
                reverse_primers.append(r_primer)

            primer_sets.append(PrimerSet(forward_primers, reverse_primers))

        return primer_sets

    def get_hetero_seqs(self, incomplete_forward_primer: MBPrimerBuilder,
                        incomplete_reverse_primer: MBPrimerBuilder,
                        forward_spacer: SpacerAlignment,
                        reverse_spacer: SpacerAlignment,
                        num_to_return: int) -> List[PrimerSet]:
        """Generates <num_to_return> PrimerSets, trying to minimise
        complementarity between the heterogeneity regions in the primers and
        that of other sequences."""
        # Generate semi-random heterogeneity seqs for the for and rev primers.
        forward_spacer_seqs = gen_hetero_set(incomplete_forward_primer,
                                             forward_spacer,
                                             self._random_per_align)
        reverse_spacer_seqs = gen_hetero_set(incomplete_reverse_primer,
                                             reverse_spacer,
                                             self._random_per_align)

        # Find the best combinations of forward and reverse primers with the
        # least spacer region complementarity.
        return self._spacer_sorter.filter_and_make_primer_sets(
            incomplete_forward_primer,
            incomplete_reverse_primer,
            forward_spacer_seqs,
            reverse_spacer_seqs,
            num_to_return)

    def set_pairwise(self, degen: bool = None) -> None:
        """Sets the primer types to be returned to pairwise primers"""
        self._spacer_sorter = SortForPairwise(self._max_spacer_length,
        self._num_hetero, num_pairings_to_comp=self._num_pairings_to_compare,
                                              degen=degen)

    def set_rigour(self, rigour: int) -> None:
        """Increases the sample sizes of the random generation process
        according to <rigour>. """
        self._rigour = rigour
        if rigour > 0:
            rigour += 1
            self._random_per_align = INITIAL_PRIMER_SET_SIZE * rigour
            # Square root to slow rate of growth.
            self._num_pairings_to_compare = int(NUM_PAIRINGS_TO_COMP *
                                             rigour ** (1/3))
            self._spacer_sorter.set_num_pairings_to_compare(
                NUM_PAIRINGS_TO_COMP * rigour)
        elif rigour < 0:
            rigour -= 1
            self._random_per_align = int(INITIAL_PRIMER_SET_SIZE /
                                         -rigour)
            # Square root to slow rate of growth.
            self._num_pairings_to_compare = int(NUM_PAIRINGS_TO_COMP /
                                             (-rigour) ** (1/3))
            self._spacer_sorter.set_num_pairings_to_compare(
                self._num_pairings_to_compare)
        elif rigour == 0:
            self._random_per_align = INITIAL_PRIMER_SET_SIZE
            # Square root to slow rate of growth.
            self._num_pairings_to_compare = NUM_PAIRINGS_TO_COMP
            self._spacer_sorter.set_num_pairings_to_compare(
                self._num_pairings_to_compare)

