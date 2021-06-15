from typing import List, Tuple
from Bio.Seq import Seq
from primer_tools import MBPrimer, MBPrimerBuilder
from presenters import Presenter, ConsolePresenter
import random

NUM_SPACERS = 4
# The default maximum length of any spacer created
MAX_SPACER_LENGTH = 12
# The default size of the heterogeneity region to proceed the binding region of
# the primer
NUM_HETERO = 12


class HeteroGen:
    """Manages the the creation of heterogeneity spacers.

    === Private Attributes ===
    _max_spacer_length:
            The maximum length of a spacer produced.
    _num_hetero:
            The length of the heterogeneity that should be ensured across."""

    _max_spacer_length: int
    _num_hetero: int
    _presenter: Presenter
    _num_spacers: int

    def __init__(self, max_spacer_length: int = MAX_SPACER_LENGTH, num_hetero: int = 12,
                 presenter: Presenter = ConsolePresenter()) -> None:
        """Initialises the attributes to the values specified."""
        self._presenter = presenter

        self._max_spacer_length = max_spacer_length
        self._num_hetero = num_hetero
        self._num_spacers = NUM_SPACERS

    def _generate_hetero_seq(self, primer: Seq, vary_len: int,
                             num_sets: int) -> Tuple[Seq]:
        """Generates heterogeneity sequences given a <primer> sequence.
        Ensures nucleotide diversity among the first <vary_len> bases in the
        sequence. Returns a list of <num_sets> tuples, each containing 4
        tandem heterogeneity sequences. Precondition: len(primer) > vary_len.
        """

        pass

    # Begin methods related to generating spacer region alignments

    def get_all_spacer_combos(self, seq: Seq) -> List[Tuple[int, int, int, int]]:
        """Given a <seq>, will provide a set of alignments of that sequence (
        produced by shifting it to the right) that ensures nucleotide
        diversity across the first <self.num_hetero> bases. Returns a list of
        tuples, each of which contain unique combinations of valid heterogeneity
        spacer lengths."""

        valid_spacer_combos = []
        for i in range(0, self._num_hetero):
            # Seeding with a single spacer length
            spacers = [i]
            self._get_all_compatible_spacers(seq, spacers, self._num_spacers,
                                             spacer_combo_list=valid_spacer_combos)
        return valid_spacer_combos

    def _get_all_compatible_spacers(self, seq: Seq, spacers: List[int],
                                    target_depth: int,
                                    depth: int = 0,
                                    spacer_combo_list: List[
                                        Tuple[int, int, int, int]] = None) -> None:
        """Produces a tree of valid combinations of spacer lengths, returns
        the first node in this tree. If a depth of <target_depth> cannot be
        reached, returns a node with node.val -1. If <spacer_combo_list> is
        specified then all valid combinations of spacers will be appended to
        it. """

        if depth == target_depth - 1:
            if spacer_combo_list is not None:
                spacer_combo_list.append(tuple(spacers))
            return None

        for spacer in self._get_compatible_spacers(seq, spacers):
            # The set of spacers to be tested for heterogeneity across
            # self.num_hetero
            trial_spacers = spacers.copy()
            trial_spacers.append(spacer)
            self._get_all_compatible_spacers(seq, trial_spacers, target_depth,
                                             depth + 1, spacer_combo_list)
        return None

    def _get_compatible_spacers(self, seq: Seq,
                                seqs_spacers: List[int]) -> Tuple[int, ...]:
        """Returns a tuple containing all spacers < <self.num_hetero> such
        that for all j + spacer < <self.num_hetero>, seq[j + spacer] != any
        seqs[i][j]. Returns an empty tuple if no such spacers exist. TODO

        """

        valid_spacers = []
        # Assume that there exist no valid spacers less than the greatest
        # already used
        start_point = seqs_spacers[-1] + 1

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
                if seq[i + spacer - spacers[j]] == seq[i]:
                    return False

        return True

    def visualise_seq_arr(self, seq_arrs: List[List[List[str or Seq]]]) -> None:
        """Displays a visual representation of a given <seq_arrs>.

        Precondition: For all valid i, j, and k,
        len(<seq_arr>[i][j][k]) > self._num_hetero"""
        to_print = ''

        for i in range(len(seq_arrs)):
            to_print += ''.join(["Sequence array #", str(i), '\n'])

            for j in range(1, min(self._num_hetero, 11)):
                to_print += str(j) + ' '
            to_print += '\n'

            # Add spaces between the bases in the heterogeneity region
            for j in range(self._num_spacers):
                to_print += ' '.join(seq_arrs[i][j][0:self._num_hetero]) + ' '
                to_print += ''.join(seq_arrs[i][j][self._num_hetero:]) + '\n'

            to_print += '\n'

        self._presenter.print(to_print)


    def visualise_spacer_combos(self, spacers: List[Tuple[int, int, int, int]],
                                seq: Seq) -> None:
        """Displays a visual representation of all of the <spacers> for the
        given <seq>. Assumes a set of 4 spacers per tuple."""

        to_print = ''
        for i in range(0, len(spacers)):
            spacer_tup = spacers[i]
            to_print += 'Spacer #' + str(i) + ' ' + str(spacer_tup) + '\n'

            for j in range(1, min(self._num_hetero, 11)):
                to_print += str(j) + ' '
            to_print += '\n'

            for j in range(0, self._num_spacers):
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

    """ Begin methods that create sequences to occupy a heterogeneity alignment.
     Main idea here is to generate a several possible combinations of spacers
     and then evaluate them based on a set of parameters."""

    def gen_hetero_seqs(self, incomplete_forward_primer: MBPrimerBuilder,
                        incomplete_reverse_primer: MBPrimerBuilder,
                        spacers: Tuple[int, int, int, int], max_GC: int) -> \
            List[Tuple[MBPrimer]]:

        pass

    def _gen_sequence_array(self, binding_seq: Seq,
                            spacers: Tuple[int, int, int, int]) -> List[List[str]]:
        """Returns a representation of the <spacers> and <binding_seq>.

        >>> HeteroGen._gen_sequence_array(Seq('ATCG'), (0, 1, 2, 3))
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
            for j in range(0, self._num_hetero - spacers[i] + 1):
                if j < len(binding_seq):
                    sequence_array[i].append(binding_seq[j])
                else:
                    sequence_array[i].append('-')

        return sequence_array

    def get_vacant_bases(self, spacers: Tuple[int, int, int, int]) -> List[
        List[int]]:
        """Returns a List containing a List of indices of unfilled bases in
        <sequence_array>.
        >>> seq_arr = [['A', 'T', 'C', 'G'],
        ...             ['',  'A','T', 'C'],
        ...             ['',  '', 'A', 'T'],
        ...             ['',  '',  '', 'A']]
        >>> HeteroGen.get_vacant_bases(seq_arr)
         [[1, 2, 3], [2, 3], [3], []]"""
        unfilled_bases = []
        for i in range(self._num_hetero + 1):
            unfilled_bases.append([])

        for i in range(self._num_spacers):
            for j in range(spacers[i]):
                unfilled_bases[j].append(i)

        return unfilled_bases

    def get_potential_bases(self, sequence_array: List[List[str]],
                            column: int) -> List[str]:
        """Returns  a list of all the bases that, if inserted into column, would maintain
        heterogeneity along the aligned nucleotide column.
        >>> seq_arr = [['A', 'T', 'C', 'G'],
        ...             ['',  'A','T', 'C'],
        ...             ['',  '', 'A', 'T'],
        ...             ['',  '',  '', 'A']]
        >>> HeteroGen.get_potential_bases(seq_arr, 0)
         ['T', 'C', 'G']
         >>> HeteroGen.get_potential_bases(seq_arr, 2)
         ['G']"""
        possible_bases = ['A', 'T', 'C', 'G']
        for i in range(self._num_spacers):
            if sequence_array[i][column]:
                possible_bases.remove(sequence_array[i][column])
        return possible_bases

    def select_and_set(self, potential_bases: List[str],
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

    def gen_heterogeneity_spacers_rand(self, binding_seq: Seq,
                                  spacers: Tuple[int, int, int, int]) -> List[Seq]:
        """Returns a list of randomly generated valid heterogeneity spacers.
        When possible, avoids base repeats on same strand."""

        sequence_array = self._gen_sequence_array(binding_seq, spacers)
        unfilled_bases = self.get_vacant_bases(spacers)

        # Iterate through each column in reverse.
        for column in range(len(unfilled_bases) - 2, -1, -1):
            potential_bases = self.get_potential_bases(sequence_array, column)
            # Randomly fill one of the unspecified bases
            while unfilled_bases[column]:
                self.select_and_set(potential_bases, unfilled_bases, column,
                                    sequence_array)

        # Extract newly generated heterogeneity spacers from the sequence array.
        final_spacers = []
        for i in range(self._num_spacers):
            spacer = Seq(str.join('', sequence_array[i][0:spacers[i]]))
            final_spacers.append(spacer)

        return final_spacers






