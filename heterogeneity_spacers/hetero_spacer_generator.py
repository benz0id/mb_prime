from typing import List, Tuple, Dict
from Bio.Seq import Seq

from presenters.presenters import Presenter, ConsolePresenter


class Node:
    """A simple class used in constructing trees of integers.
    === Public Attributes ===
    val:
        The integer value contained by this node.
    nexts:
        A list of nodes pointed to.
    """

    val: int

    def __init__(self) -> None:
        """Sets attributes to default values"""
        self.val = -1
        self.next_list = []
        self.parent = None


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
                 presenter: Presenter = None) -> None:
        """Initialises the attributes to the values specified."""
        if presenter is None:
            self._presenter = ConsolePresenter()

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
            self._build_spacer_tree(seq, spacers, 4,
                                    spacer_combo_list=valid_spacer_combos)
        return valid_spacer_combos

    def _build_spacer_tree(self, seq: Seq, spacers: List[int],
                           target_depth: int,
                           depth: int = 0,
                           spacer_combo_list: List[Tuple[int]] = None) -> Node:
        """Produces a tree of valid combinations of spacer lengths, returns
        the first node in this tree. If a depth of <target_depth> cannot be
        reached, returns a node with node.val -1. If <spacer_combo_list> is
        specified then all valid combinations of spacers will be appended to
        it. """

        # TODO remove nodes

        node = Node()

        if depth == target_depth - 1:
            node.val = spacers[-1]
            if spacer_combo_list is not None:
                spacer_combo_list.append(tuple(spacers))
            return node

        for spacer in self._get_compatible_spacers(seq, spacers):
            # The set of spacers to be tested for heterogeneity across
            # self.num_hetero
            trial_spacers = spacers.copy()
            trial_spacers.append(spacer)
            nxt = self._build_spacer_tree(seq, trial_spacers,
                                          target_depth,
                                          depth + 1, spacer_combo_list)
            if nxt.val > 0:
                node.next_list.append(nxt)

        if len(node.next_list) > 0:
            node.val = spacers[-1]
        return node

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
                ind = i + spacer - spacers[j]  # TODO REMOVE
                if seq[ind] == seq[i]:
                    return False

        return True

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
