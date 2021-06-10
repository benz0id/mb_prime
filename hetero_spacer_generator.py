from typing import List, Tuple
from Bio.Seq import Seq
from primer_tools import MBPrimer, MBPrimerBuilder
from presenters import Presenter, ConsolePresenter
import random


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

    def __init__(self, max_spacer_length: int = 12, num_hetero: int = 12,
                 presenter: Presenter = ConsolePresenter()) -> None:
        """Initialises the attributes to the values specified."""
        self._presenter = presenter

        self._max_spacer_length = max_spacer_length
        self._num_hetero = num_hetero

    def _generate_hetero_seq(self, primer: Seq, vary_len: int,
                             num_sets: int) -> Tuple[Seq]:
        """Generates heterogeneity sequences given a <primer> sequence.
        Ensures nucleotide diversity among the first <vary_len> bases in the
        sequence. Returns a list of <num_sets> tuples, each containing 4
        tandem heterogeneity sequences. Precondition: len(primer) > vary_len.
        """

        pass

    # Begin methods related to generating spacer region alignments

    def get_all_spacer_combos(self, seq: Seq) -> List[Tuple[int]]:
        """Given a <seq>, will provide a set of alignments of that sequence (
        produced by shifting it to the right) that ensures nucleotide
        diversity across the first <self.num_hetero> bases. Returns a list of
        tuples, each of which contain unique combinations of valid heterogeneity
        spacer lengths."""

        valid_spacer_combos = []
        for i in range(0, self._num_hetero):
            # Seeding with a single spacer length
            spacers = [i]
            self._get_all_compatible_spacers(seq, spacers, 4,
                                             spacer_combo_list=valid_spacer_combos)
        return valid_spacer_combos

    def _get_all_compatible_spacers(self, seq: Seq, spacers: List[int],
                                    target_depth: int,
                                    depth: int = 0,
                                    spacer_combo_list: List[
                                        Tuple[int]] = None) -> None:
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
                                seqs_spacers: List[int]) -> Tuple[int]:
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

    def visualise_spacer_combos(self, spacers: List[Tuple[int]],
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

            for j in range(0, 4):
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
                        spacers: Tuple[int], max_GC: int) -> \
            List[Tuple[MBPrimer]]:

        pass

    def _gen_sequence_array(self, binding_seq: Seq,
                            spacers: Tuple[int]) -> List[List[str]]:
        """Returns a representation of the <spacers> and <binding_seq>.

        >>> HeteroGen._gen_sequence_array(Seq('ATCG'), (0, 1, 2, 3))
            [['A', 'T', 'C', 'G'],
            ['',  'A','T', 'C'],
            ['',  '', 'A', 'T'],
            ['',  '',  '', 'A']]"""
        sequence_array = []  # [row][column]
        for i in range(len(spacers)):
            sequence_array[i] = []
            # Fill hetero spacer region with empty slots
            for _ in range(spacers[i]):
                sequence_array[i].append('')  # Indicates an empty slot

            # Occupy remaining region (plus an extra base) with primer
            # Third spacer is always the largest.
            for j in range(0, len(sequence_array[i]) - spacers[3] + 1):
                sequence_array[i].append(binding_seq[j])
        return sequence_array

    def get_vacant_bases(self, sequence_array: List[List[str]]) -> List[
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
        for i in range(len(sequence_array[1])):
            unfilled_bases.append(0)
            for j in range(4):
                unfilled_bases[i] = []
                if sequence_array[j][i] == '':
                    unfilled_bases[i].append(j)
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
        for i in range(4):
            if sequence_array[column][i]:
                possible_bases.remove(sequence_array[column][i])
        return possible_bases

    def gen_heterogeneity_spacers(self, binding_seq: Seq,
                                  spacers: Tuple[int]) -> List[Seq]:

        sequence_array = self._gen_sequence_array(binding_seq, spacers)
        unfilled_bases = self.get_vacant_bases(sequence_array)

        # Iterate through each column in reverse.
        for i in range(len(unfilled_bases), -1, -1):
            potential_bases = self.get_potential_bases(sequence_array, i)
            # Randomly fill one of the unspecified bases
            while unfilled_bases[i]:
                specific_potential_bases = potential_bases.copy()
                base_index = random.choice(unfilled_bases[i])

                # Avoid repeating the last base in the sequence if possible
                previous_base = sequence_array[i + 1][base_index]
                if len(specific_potential_bases) > 1 and \
                        previous_base in specific_potential_bases:
                    specific_potential_bases.remove(previous_base)

                selected_base = random.choice(specific_potential_bases)

