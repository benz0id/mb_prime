import functools
from abc import ABC, abstractmethod
from typing import List, Tuple, TypeVar, Union
from Bio.Seq import Seq
from multiprocessing import Process
from bisect import bisect_right
from hetero_spacer_generator.primer_tools import Comparable, get_n_highest_sbs
from Bio.Align import MultipleSeqAlignment
from Bio.Align.AlignInfo import SummaryInfo

# Where to change if new parameters are found:
#   Binding Params
#
#
#

C = TypeVar('C')
DO_ERROR_CHECKING = True


def reverse_inds(inds: List[int], size: int,
                 err_chk: bool = DO_ERROR_CHECKING) -> List[int]:
    """If <ind> is some number of indices on a list of <size>, this function
    returns the location of those indices if the list were to be reversed.

    Precondition:
        0 <= each <ind> <= size"""

    rev_inds = []
    if err_chk:
        err_str = ''
        if size < 0:
            err_str = "Received negative hypothetical array size"
        for ind in inds:
            if size < ind:
                err_str = "Received index outside of hypothetical array size."
                break
        if err_str:
            raise ValueError(err_str)

    for ind in inds:
        rev_inds.append(size - ind - 1)
    return rev_inds


@functools.total_ordering
class Param(Comparable):
    """Stores some parameters and their associated overall score.

    _score:
        The score of this parameter"""

    _score: float

    def __init__(self) -> None:
        self._score = 0

    def __lt__(self, other) -> bool:
        return self._score < other.get_score()

    def __eq__(self, other) -> bool:
        return self._score == other.get_score()

    def set_score(self, score: float) -> None:
        """Sets this params score to the given value."""
        self._score = score

    def get_score(self) -> None:
        """Returns this param's score."""
        print(self._score)


class BindingParams(Param):
    """The parameters used to construct a single binding sequence.

    === Public Attributes ===

    _5p_ind:
        The index of the 5' base of this primer.
    _len:
        The length of this sequence.
    _adapter_num:
        The adapter sequence with which this sequence is paired."""

    _5p_ind: int
    _len: int
    _adapter_num: int

    def __init__(self, fivep_ind: int, length: int, adapter_num: int = -1) \
            -> None:
        """Initialises the this parameters object with the given values."""
        super().__init__()
        self._5p_ind = fivep_ind
        self._len = length
        self._adapter_num = adapter_num

    def get_5p_ind(self) -> None:
        """Returns this param's 5p index."""
        print(self._5p_ind)

    def get_len(self) -> None:
        """Returns this param's length."""
        print(self._len)

    def set_adapter_num(self, adapter_num: int) -> None:
        """Sets this param's adapter number."""
        self._adapter_num = adapter_num

    def get_adapter_num(self) -> None:
        """Returns this param's adapter number."""
        print(self._adapter_num)


class BindingPairParams(Param):
    """The parameters used to construct a pair of forward and reverse binding
    sequences.

    === Public Attributes ===
    f_params:
        Parameters of the forward binding seq.
    r_params:
        Parameters of the reverse binding seq."""

    f_params: BindingParams
    r_params: BindingParams

    def __init__(self, f_params: BindingParams, r_params: BindingParams) \
            -> None:
        """Initialises this class with the given foward and reverse binding
        parameters."""
        super().__init__()
        self.f_params = f_params
        self.r_params = r_params


P = TypeVar('P', Param)


class BindingSeqIterator(ABC):
    """A class that iterates over all possible permutations of a given set of
    parameters of binding sequences, returning the binding sequences specified
    by those parameters.

    === Private Attributes ===
    _consensus:
        The sequence from which binding sequences will be extracted.
    """

    @abstractmethod
    def __next__(self):
        pass

    def __iter__(self):
        return self

    @abstractmethod
    def get_last_param(self) -> Param:
        """Returns the parameters of the binding sequence returned by the last
        call as a Param object"""
        pass


class BindingSeqScorer(ABC):
    """A class designed to score and collect the best outputs of a
    BindingSeqIterator.

    === Private Attributes ===
    _best_params:
        A list containing the highest scoring parameters. Sorted L-G,
        Best-Worst.
    _seq_iterator:
        An iterator that produces some set of binding sequences, capable of
        providing the generation parameters for those binding sequences if
        requested.
    _min_score:
        The least score of a parameter in _best_params.
    _num_params_to_keep:
        The length of _best_params
    _cur_num_params:
        The number of non-empty indices in _best_params
    _consensus:
        The consensus sequence from which this class will extract binding
        Sequences.
        """

    _best_params: List[P]
    _seq_iterator: BindingSeqIterator
    _min_score: float
    _num_params_to_keep: int
    _cur_num_params: int
    _consensus: Seq

    def __init__(self, num_params_to_keep: int, consensus: Seq) -> None:
        """Initialises this scorer with the appropriate default values."""
        self._best_params = []
        self._min_score = 0
        self._cur_num_params = 0
        self._num_params_to_keep = num_params_to_keep
        self._consensus = consensus

    def _add_to_best_params(self, param: P) -> None:
        """Adds <param> to <self._best_params>, maintaining G-L order."""
        ind = bisect_right(self._best_params, param)
        self._best_params.insert(ind, param)
        # Remove last ind if desired length has been reached.
        if len(self._best_params) > self._num_params_to_keep:
            self._best_params.pop(-1)
        else:
            self._cur_num_params += 1

    def _add_if_lower(self, score: float) -> None:
        """Adds the parameters used to construct the last binding sequence(s) to
        self._best_params iff it is better than the worst score in the list or
        the list is incomplete."""
        if score < self._min_score or \
                self._cur_num_params < self._num_params_to_keep:
            param = self._seq_iterator.get_last_param()
            self._add_to_best_params(param)

    @classmethod
    def get_best_params(cls) -> List[P]:
        """Returns the best params found by this class."""
        return cls._best_params

    def _split_into_n(self, num_children: int) -> List[P]:
        """Finds and returns the parameters that produce the best binding
        sequences by splitting the task among <num_threads> children."""
        split_classes = self._split_class(num_children)
        children = []
        best_params_comp = []
        for bss in split_classes:
            children.append(Process(target=bss.find_best))
        for child in children:
            child.start()
        for child in children:
            child.join()
        for splt in split_classes:
            best_params_comp.append(splt.get_best_params())
        self._best_params = get_n_highest_sbs(best_params_comp,
                                              self._num_params_to_keep,
                                              least=True)
        return self._best_params

    @abstractmethod
    def run_get_best_params(self, num_children: int = 1) -> List[P]:
        """Return the parameters that produce the best binding sequences."""
        pass

    @abstractmethod
    def _get_iterator(self) -> BindingSeqIterator:
        """Return the iterator from which this BindingSeqScorer will extract
        binding sequences."""
        pass

    @classmethod
    @abstractmethod
    def _split_class(cls: C, num: int) -> List[C]:
        """Splits the binding sequence space into <num> other BindingSeqScorers
        and returns them."""
        pass

    @abstractmethod
    def find_best(self) -> None:
        """Completely iterates through <self._seq_iterator> and stores the best
        scoring params in <self._best_params>"""
        pass


class HomoSeqIterator(BindingSeqIterator):
    """A class designed to iterate over all possible binding sequences given
    some parameters.
    _allowed_5p:
        List of allowable 5' starting indices for the binding sequence
    _allowed_lengths:
        List of all allowed lengths for the primer.
    _consensus:
        The consensus sequence from which this class will extract binding
        Sequences. 5' - 3'
    _i:
        The current iteration number. Corresponds to parameters
        _allowed_5p[self._i // len(self._allowed_5p)] and
        _allowed_5p[self._i % len(self._allowed_5p)]
    _max_iter:
        The maximum allowed value for an iteration.
    """
    _allowed_5p: List[int]
    _allowed_lengths: List[int]
    _consensus: str
    _i: int
    _max_iter: int

    def __init__(self, consensus: str or Seq, allowed_5p: List[int],
                 allowed_lengths: List[int]) -> None:
        """Initialises this class using the given values."""
        self._i = -1
        self._max_iter = len(allowed_lengths) * len(allowed_5p)
        self._allowed_5p = allowed_5p
        self._allowed_lengths = allowed_lengths
        self._consensus = str(consensus)
        if max(allowed_5p + [0]) + max(allowed_lengths + [0]) > len(consensus):
            raise ValueError("One or more specified primers exceed the bounds "
                             "of the given consensus")

    def __next__(self) -> str:
        """Gets the next binding sequence, 5' - 3'"""

        self._i += 1
        if self._max_iter == self._i:
            raise StopIteration

        cur_5p = self._allowed_5p[self._i // len(self._allowed_lengths)]
        cur_len = self._allowed_lengths[self._i % len(self._allowed_lengths)]

        # Return specified portion of alignment.
        return self._consensus[cur_5p: cur_5p + cur_len]

    def get_last_param(self) -> BindingParams:
        """Returns the binding parameters of the sequence last returned."""
        last_5p = self._allowed_5p[self._i // len(self._allowed_5p)]
        last_len = self._allowed_lengths[self._i % len(self._allowed_5p)]
        return BindingParams(last_5p, last_len)

    def get_last_len(self) -> int:
        """Returns the length of the last returned binding sequence."""
        return self._allowed_lengths[self._i % len(self._allowed_5p)]

    def get_last_5p(self) -> int:
        """Returns the length of the last returned binding sequence."""
        return self._allowed_5p[self._i // len(self._allowed_lengths)]


class HomoSeqScorer(BindingSeqScorer):
    """Evaluates the features of primers produced by all possible
    permutations of some construction parameters.

    for_binding |       5'-GCATGGTGATCGT-3'
    <consensus> | 5'-ATCGATGCATGGTGATCGTAGGTCGGTAGTCGTCAGTCGTGTAGTCGG-3'
    rev_binding |                              3'-AGCAGTCAGCACA-5'

    === Private Attributes ===

    _alignment:
        The alignment for which binding sequences will be analysed.
    _allowed_5p:
        List of allowable 5' starting indices for the binding sequence
    _allowed_lengths:
        List of all allowed lengths for the primer.
    _adapters:
        A list of adapter sequences 5' - 3'.
    _consensus:
        The consensus sequence of the alignment.
    _rev:
        Whether this scorer is considering reverse primers.
    """
    _alignment: MultipleSeqAlignment
    _allowed_5p: List[int]
    _allowed_lengths: List[int]
    _adapters: List[str]
    _consensus: str
    _rev: bool

    def __init__(self, alignment: MultipleSeqAlignment, allowed_5p: List[int],
                 allowed_lengths: List[int], adapters: List[str],
                 num_params_to_keep: int, rev: bool) -> None:
        """Constructs """
        # TODO Find smarter ways of extracting consensus.
        sa = SummaryInfo(alignment)
        consensus = sa.dumb_consensus()
        # If considering reverse primers, work on reverse complement of
        # consensus.
        if rev:
            consensus = consensus.reverse_complement()
            # Accounting for reversed positions.
            self._allowed_5p = reverse_inds(allowed_5p, len(consensus))
        else:
            self._allowed_5p = allowed_5p

        super().__init__(num_params_to_keep, consensus)

        self._consensus = str(consensus)
        self._alignment = alignment
        self._allowed_lengths = allowed_lengths
        self._adapters = adapters
        self._rev = rev

    def _get_iterator(self) -> BindingSeqIterator:
        """Returns an iterator over all possible binding sequences."""
        return HomoSeqIterator(consensus=self._consensus,
                               allowed_5p=self._allowed_5p,
                               allowed_lengths=self._allowed_lengths)

    @classmethod
    def _split_class(cls: C, num: int) -> List[C]:
        """Splits the class into <num> separate homo_seq scorers, each sampling
        the same number of binding sequences. To be used in threading."""
        pass

    def score_binding_seq(self, adapter_ind: int, binding_seq: str) -> float:
        """Returns the score of the given adapter-binding seq combo."""
        pass

    def find_best(self) -> None:
        """Finds the binding sequences with the lowest scores and stores them
        in <self._best_params>"""
        pass

    def run_get_best_params(self, num_children: int = 1) -> List[BindingParams]:
        """Iterates over all possible binding sequences and returns the best
        results. Iff num_children > 1, spawns <num_children> processes and
        splits the binding space between them."""
        if num_children > 1:
            return self._split_into_n(num_children)


class HeteroSeqIterator(BindingSeqIterator):
    """A class designed to iterate over all possible pairs of binding sequences
    given some parameters

    _consensus:
        The sequence from which potential pairs of primers will be extracted.
    _rc_consensus:
        The reverse compliment of the consensus, used to generate reverse
        primers.

    _f_allowed_5p:
        The list of allowable 5p or starting indices of any given forward
        binding sequence (relative to original consensus).
    _r_allowed_5p:
        The list of allowable 5p or starting indices of any given reverse primer
        (relative to consensus, to reverse compliment consensus).

    _f_allowed_len:
        The list of allowable lengths for any returned forward binding sequence.
    _r_allowed_len:
        The list of allowable lengths for any returned reverse binding sequence.

    _allowed_amp_len:
        The list of allowable amplicon lengths, forward 5' to reverse 5'.
    _min_allowed_amp_len:
        The least allowed amplicon length

    _for_iter:
        Responsible for iterating over every possible forward binding sequence.
    _rev_iter:
        Responsible for iterating over a subset of allowable reverse binding
        sequences. Subset is decided by current state of for_iter.
    _cur_for:
        The current forward binding sequence.
    """
    _consensus: str
    _rc_consensus: str

    _f_allowed_5p: List[int]
    _r_allowed_5p: List[int]

    _f_allowed_len: List[int]
    _r_allowed_len: List[int]

    _allowed_amp_len: List[int]
    _min_allowed_amp_len: int

    _for_iter: HomoSeqIterator
    _rev_iter: HomoSeqIterator
    _cur_for: str

    def __init__(self, consensus: Union[Seq, str], f_allowed_5p: List[int],
                 r_allowed_5p: List[int], f_allowed_len: List[int],
                 r_allowed_len: List[int], allowed_amp_len: List[int]) -> None:
        """Initialises a HeteroSeqIterator that iterates over all primers that
        match the above criteria. All 5p indices are relative to consensus,
        not its reverse compliment.

        Preconditions:
            Binding Seqs cannot be longer than amplicon:
            min(allowed_amp_length) > max(r_allowed_len) + max(f_allowed_len)
            """
        consensus_cpy = Seq(str(consensus[:]))
        self._rc_consensus = str(consensus_cpy.reverse_complement())
        self._consensus = str(consensus_cpy)

        self._f_allowed_5p = f_allowed_5p
        self._r_allowed_5p = r_allowed_5p

        self._f_allowed_len = f_allowed_len
        self._r_allowed_len = r_allowed_len

        self._allowed_amp_len = allowed_amp_len
        self._min_allowed_amp_len = min(allowed_amp_len)

        # Configure attributes in preparation for first iteration.
        self._for_iter = HomoSeqIterator(self._consensus, f_allowed_5p,
                                         f_allowed_len)
        self._update_homo_iterators()

    def __next__(self) -> Tuple[str, str]:
        """Gets the next potential forward reverse primer binding seqs formatted
        as strings, both 5' - 3'."""

        # Try to fetch the next potential reverse binding seq.
        try:
            cur_rev = self._rev_iter.__next__()
        except StopIteration:
            # We've run out of potential reverse binding seqs for the current
            # forward binding seq.
            cur_rev = ''

        # If we've run out of Increment forward binding seq
        while cur_rev == '':
            self._update_homo_iterators()
            if self._cur_for == '':
                # We're out of potential forward binding seqs.
                raise StopIteration
            try:
                cur_rev = self._rev_iter.__next__()
            except StopIteration:
                pass

        return self._cur_for, cur_rev

    def _update_homo_iterators(self) -> None:
        """To be called when the reverse iterator runs out of valid binding
        sequences for the current forward binding sequence.

        Increments the forward HomoSeqGenerator and creates a new reverse
        HomoSeqIterator that contains all valid reverse binding sequences for
        the new forward binding sequence."""
        try:
            self._cur_for = self._for_iter.__next__()
        except StopIteration:
            # We've run out of valid forward binding sequences.
            self._cur_for = ''
            return

        f_5p = self._for_iter.get_last_5p()

        allowed_rev_5p = self._get_allowed_rev_5p(f_5p)
        allowed_rev_5p_rc = reverse_inds(allowed_rev_5p, len(self._consensus))
        self._rev_iter = HomoSeqIterator(self._rc_consensus, allowed_rev_5p_rc,
                                         self._r_allowed_len)

    def _get_allowed_rev_5p(self, f_5p: int) -> List[int]:
        """Returns a list of all allowable reverse 5p inds relative to consensus
        sequence."""
        allowed = []
        # Which 5' binding sites can produce an amplicon of a desired size?
        for r_5p in self._r_allowed_5p:
            # If the reverse 5' site is too small to form any valid amplicon
            # continue.
            if r_5p + 1 < f_5p + self._min_allowed_amp_len:
                continue

            for amp_len in self._allowed_amp_len:
                if r_5p + 1 == amp_len + f_5p:
                    # This 5' binding site produces an amplicon of the desired
                    # length. We don't need to check any others so break.
                    allowed.append(r_5p)
                    break
        return allowed

    def get_last_param(self) -> BindingPairParams:
        """Returns the binding parameters of the sequence last returned."""
        return BindingPairParams(self._for_iter.get_last_param(),
                                 self._rev_iter.get_last_param())


class HeteroSeqScorer(BindingSeqScorer):
    """
    Evaluates the features of each allowable pair of forward and reverse
    binding sequence.

    _consensus:

    _f_allowed_5p:
    _r_allowed_5p:

    _f_allowed_len:
    _r_allowed_len:

    _allowed_amp_len:

    _adapters

 """

    def __init__(self, consensus: str, f_allowed_5p: List[int],
                 r_allowed_5p: List[int], f_allowed_len: List[int],
                 r_allowed_len: List[int], allowed_amp_len: List[int],
                 adapters: List[str]) -> None:
        """Initialises a hetero"""

    def run_get_best_params(self, num_children: int = 1) -> List[Param]:
        pass

    def _get_iterator(self) -> BindingSeqIterator:
        pass

    @classmethod
    def _split_class(cls: C, num: int) -> List[C]:
        pass

    def find_best(self) -> None:
        pass


class BestPrimers:
    """A class intended to rigourously search for the best possible primers in
    a given region of an alignment.

    === Private Attributes ===

    _alignment:
        The alignment for which binding sequences will be analysed.

    _f_region_start:
        The starting index of the
    _r_region_start:
    _f_region_end:
    _r_region_end:

    _f_adapters:
    _r_adapters:

    _f_max_length:
    _r_max_length:
    _f_min_length:
    _r_min_length:

    _min_amp_len:
    _max_amp_len:


    """
