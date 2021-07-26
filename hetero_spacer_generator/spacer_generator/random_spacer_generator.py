import random
from abc import abstractmethod
from typing import Dict, List, Tuple
from ..primer_types import SpacerSet
from Bio.Seq import Seq

from hetero_spacer_generator.defaults import INITIAL_PRIMER_SET_SIZE, \
    NUM_PAIRINGS_TO_COMP, NUM_HETERO, NUM_SPACERS
from hetero_spacer_generator.spacer_generator.spacer_filters import \
    SortForSimultaneous
from hetero_spacer_generator.primer_tools import HeteroSeqTool, \
    IncompatibleSpacerError, MBPrimerBuilder, \
    PrimerSet
from hetero_spacer_generator.primer_types import SpacerAlignment


def get_vacant_bases(spacers: SpacerAlignment) -> List[List[int]]:
    """Returns a List containing a List of indices of unfilled bases in
    <sequence_array>."""
    unfilled_bases = []
    for i in range(spacers[3] + 1):
        unfilled_bases.append([])

    for i in range(NUM_SPACERS):
        for j in range(spacers[i]):
            unfilled_bases[j].append(i)

    return unfilled_bases


def get_potential_bases(sequence_array: List[List[str]],
                        column: int) -> List[str]:
    """Returns  a list of all the bases that, if inserted into column, would
    maintain heterogeneity along the aligned nucleotide column."""
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


def gen_heterogeneity_spacers_rand(binding_seq: Seq,
                                   spacers: SpacerAlignment):
    """Returns a list of randomly generated valid heterogeneity spacers.
    When possible, avoids base repeats on same strand."""

    sequence_array = gen_sequence_array(binding_seq, spacers)
    unfilled_bases = get_vacant_bases(spacers)

    # Iterate through each column in reverse.
    for column in range(len(unfilled_bases) - 2, -1, -1):
        potential_bases = get_potential_bases(sequence_array, column)
        # Randomly fill one of the unspecified bases
        while unfilled_bases[column]:
            select_and_set(potential_bases, unfilled_bases, column,
                           sequence_array)

    # Extract newly generated heterogeneity spacers from the sequence array.
    final_spacers = []
    for i in range(NUM_SPACERS):
        spacer = Seq(str.join('', sequence_array[i][0:spacers[i]]))
        final_spacers.append(spacer)

    # Assumes NUM_SPACERS is 4
    return tuple(final_spacers)


def gen_hetero_set(incomplete_primer: MBPrimerBuilder,
                   spacers: SpacerAlignment, num_to_gen: int) \
        -> List[SpacerSet]:
    """Generates <num_to_gen> semi-random heterogeneity spaces.
    Ensures that primers match the lengths specified in <spacers> and ensure
    heterogeneity with the <incomplete_primer>'s binding_seq."""

    random_spacers = []
    binding_seq = incomplete_primer.get_binding_seq()
    for i in range(num_to_gen):
        random_spacers.append(
            gen_heterogeneity_spacers_rand(binding_seq, spacers))

    return random_spacers


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
    _number_to_cross_compare: int
    _spacer_sorter: SortForSimultaneous

    def __init__(self, max_spacer_length: int, num_hetero: int,
                 rigour: int = 1) -> None:
        """Initialises helper classes and attributes."""
        super().__init__(max_spacer_length, num_hetero)
        self._spacer_sorter = SortForSimultaneous(max_spacer_length,
                                                  num_hetero)
        self.set_rigour(rigour)

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

    def set_rigour(self, rigour: int) -> None:
        """Increases the sample sizes of the random generation process
        according to <rigour>. """
        if rigour >= 1:
            self._random_per_align = INITIAL_PRIMER_SET_SIZE * rigour
            self._number_to_cross_compare = NUM_PAIRINGS_TO_COMP * rigour
        elif rigour < 1:
            self._random_per_align = int(INITIAL_PRIMER_SET_SIZE /
                                         -rigour)
            self._number_to_cross_compare = int(NUM_PAIRINGS_TO_COMP /
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


def select_and_set(potential_bases: List[str],
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


def gen_sequence_array(binding_seq: Seq,
                       spacers: SpacerAlignment) \
        -> List[List[str]]:
    """Returns a representation of the <spacers> and <binding_seq>."""
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
