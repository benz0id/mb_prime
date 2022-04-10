import functools
from abc import ABC, abstractmethod
from typing import Dict, List, Tuple, TypeVar, Union
from Bio.Seq import Seq
from bisect import bisect_right
from hetero_spacer_generator.defaults import CONSEC_WEIGHT, CONS_WEIGHT, \
    DIMER_WEIGHT, TOTAL_WEIGHT
from hetero_spacer_generator.primer_tools import Comparable, MaxInt
from seq_alignment_analyser.align import MSA
from hetero_spacer_generator.sequence_tools import SeqAnalyzer

# Where to change if new parameters are found:
#   Binding Params
#
#
#

C = TypeVar('C')
DO_ERROR_CHECKING = True
PRINT_EVERY = 10000
seq_anal = SeqAnalyzer(degen=False)


def get_dimer_score(seq1: str, seq2: str) -> float:
    """Returns the weighted dimerisation score of the two given seqs.

    Compares sequences as though they are in a homodimer, both should be 5'-3'.

    seq1 = aaaaggA
    seq2 = ttttccT

    Demo comaprison:

        aaaaggA
        |  |  |
        Tcctttt


    """
    consec = seq_anal.comp_seqs_any_overlap(
        Seq(seq1), Seq(seq2), seq_anal.get_consec_complementarity) * \
             CONSEC_WEIGHT
    total = seq_anal.comp_seqs_any_overlap(
        Seq(seq1), Seq(seq2), seq_anal.get_non_consec_complementarity) * \
            TOTAL_WEIGHT

    return 15 * (consec + total) / (len(seq1) + len(seq2))


def get_heterodimer_score(a1: str, b1: str, a2: str, b2: str) -> float:
    """Returns the weighted dimerisation score of the  given seqs. All seqs
     should be 5'-3'."""
    score = 0
    # Score schould give a rough idea as to how stable a dimer structure would
    # be.
    score += get_dimer_score(a1, b2) ** 2
    score += get_dimer_score(b1, b2) ** 2
    score += get_dimer_score(b1, a2) ** 2
    score += get_dimer_score(a1, a2) ** 2
    return score ** (1 / 2) / 4


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

    _mean_conservation: float
    _dimer_score: float

    def __init__(self) -> None:
        pass

    def __lt__(self, other) -> bool:
        return self.get_score() < other.get_score()

    def __eq__(self, other) -> bool:
        return self.get_score() == other.get_score()

    def get_score(self) -> float:
        """Returns this param's score."""
        return ((100 - self._mean_conservation) ** CONS_WEIGHT *
                self._dimer_score ** DIMER_WEIGHT)

    def set_dimer_score(self, dimer_score: float) -> None:
        """Sets this params score to the given value."""
        self._dimer_score = dimer_score

    def get_dimer_score(self) -> float:
        """Returns this param's dimer_score."""
        return self._dimer_score

    def set_mean_conservation(self, mean_conservation: float) -> None:
        """Sets this params mean conservation to the given value."""
        self._mean_conservation = mean_conservation

    def get_mean_conservation(self) -> float:
        """Returns this param's mean conservation."""
        return self._mean_conservation


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

    def get_5p_ind(self) -> int:
        """Returns this param's 5p index."""
        return self._5p_ind

    def rev_5p_ind(self, alignment_size: int) -> None:
        """Reverses this binding parmeters position on the alignment."""
        self._5p_ind = alignment_size - self._5p_ind - 1

    def get_len(self) -> int:
        """Returns this param's length."""
        return self._len

    def set_adapter_num(self, adapter_num: int) -> None:
        """Sets this param's adapter number."""
        self._adapter_num = adapter_num

    def get_adapter_num(self) -> int:
        """Returns this param's adapter number."""
        return self._adapter_num

    def __str__(self) -> str:
        """Returns a string representation of this parameter."""
        return ''.join([
            '5\' Start site: ', str(self._5p_ind + 1), '\n',
            'Length: ', str(self._len), '\n',
            'Mean conservation: ', str(self._mean_conservation), '\n',
        ])

    def __eq__(self, other) -> bool:
        return self._5p_ind == other.get_5p_ind() and \
               self._len == other.get_len()

def get_region(param: BindingParams, rev: bool = False, dir: str = 'f') -> str:
    """Extracts the given binding param from the consensus."""
    rang = [param.get_5p_ind(),
            param.get_5p_ind() + param.get_len() - 1]

    if dir == 'r':
        rang = [param.get_5p_ind() - param.get_len() + 1, param.get_5p_ind()]

    return tuple(rang)

def binding_param_to_seq(consensus: str, param: BindingParams,
                         rev: bool = False, dir: str = 'f') -> str:
    """Extracts the given binding param from the consensus."""
    rang = [param.get_5p_ind(),
            param.get_5p_ind() + param.get_len() - 1]

    if dir == 'r':
        rang = [param.get_5p_ind() - param.get_len() + 1, param.get_5p_ind()]

    if rev:
        rang = reverse_inds(rang, len(consensus))[::-1]

    start, stop = rang
    return consensus[start: stop + 1]


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
    _final_score: int


    def __init__(self, f_params: BindingParams, r_params: BindingParams) \
            -> None:
        """Initialises this class with the given foward and reverse binding
        parameters."""
        super().__init__()
        self.f_params = f_params
        self.r_params = r_params

    def set_final_score(self, score: int) -> None:
        self._final_score = score

    def get_final_score(self) -> int:
        return self._final_score

    def __str__(self) -> str:
        return str(self.f_params) + str(self.r_params)

    def __hash__(self):
        return hash(str(self))

    def __eq__(self, other):
        return hash(other) == hash(self)


P = TypeVar('P')


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

    P = TypeVar('P')

    _best_params: List[P]
    _min_score: float
    _num_params_to_keep: int
    _cur_num_params: int
    _consensus: Seq
    _alignment: MSA

    def __init__(self, num_params_to_keep: int, consensus: Seq,
                 alignment: MSA) -> None:
        """Initialises this scorer with the appropriate default values."""
        self._best_params = []
        self._min_score = 0
        self._cur_num_params = 0
        self._num_params_to_keep = num_params_to_keep
        self._consensus = consensus
        self._alignment = alignment

    def _add_to_best_params(self, param: P) -> None:
        """Adds <param> to <self._best_params>, maintaining G-L order."""
        ind = bisect_right(self._best_params, param)
        self._best_params.insert(ind, param)
        # Remove last ind if desired length has been reached.
        if len(self._best_params) > self._num_params_to_keep:
            self._best_params.pop(-1)
        else:
            self._cur_num_params += 1

    @abstractmethod
    def get_best_params(self) -> List[Param]:
        """Returns the best params found by this class."""
        pass

    def _split_into_n(self, num_children: int):
        """Finds and returns the parameters that produce the best binding
        sequences by splitting the task among <num_threads> children."""
        """
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
        """
        raise NotImplementedError('This feature is not implemented.')

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
    def find_best(self) -> Param:
        """Completely iterates through <self._seq_iterator> and stores the best
        scoring params in <self._best_params>"""
        pass

    def _find_average_conservation(self, param: BindingParams,
                                   rev: bool, dir: str = 'f') -> None:
        """Finds and stores the average conservation across the binding region
        specified in <param>. <rev> specifies whether the primer refers to the
        sequence on the reverse complement of the consensus sequence. dir refers
        to the direction that the sequence points on in its consensus."""
        if dir == 'f':
            rang = [param.get_5p_ind(),
                    param.get_5p_ind() + param.get_len() - 1]
        else:
            rang = [param.get_5p_ind() - param.get_len() + 1,
                    param.get_5p_ind()]
        if rev:
            rang = reverse_inds(rang, len(self._consensus))[::-1]

        start, stop = rang
        cons = self._alignment.get_mean_conservation(start, stop)
        param.set_mean_conservation(cons)


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
    _cur_5p: int
    _cur_len: int

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

        self._cur_5p = self._allowed_5p[self._i // len(self._allowed_lengths)]
        self._cur_len = self._allowed_lengths[
            self._i % len(self._allowed_lengths)]

        # Return specified portion of alignment.
        return self._consensus[self._cur_5p: self._cur_5p + self._cur_len]

    def get_last_param(self) -> BindingParams:
        """Returns the binding parameters of the sequence last returned."""
        return BindingParams(self._cur_5p, self._cur_len)

    def get_last_len(self) -> int:
        """Returns the length of the last returned binding sequence."""
        return self._cur_len

    def get_last_5p(self) -> int:
        """Returns the length of the last returned binding sequence."""
        return self._cur_5p


class HomoSeqScorer(BindingSeqScorer):
    """Evaluates the features of primers produced by all possible
    permutations of some construction parameters.

    for_binding |       5'-GCATGGTGATCGT-3'
    <consensus> | 5'-ATCGATGCATGGTGATCGTAGGTCGGTAGTCGTCAGTCGTGTAGTCGG-3'

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
    _seq_iterator:
        An iterator that produces some set of binding sequences, capable of
        providing the generation parameters for those binding sequences if
        requested.
    """
    _alignment: MSA
    _allowed_5p: List[int]
    _allowed_lengths: List[int]
    _adapters: List[str]
    _seq_iterator: HomoSeqIterator
    _consensus: str
    _rev: bool

    def __init__(self, alignment: MSA, allowed_5p: List[int],
                 allowed_lengths: List[int], adapters: List[str],
                 num_params_to_keep: int, rev: bool) -> None:
        """Constructs necessary components for sampling primer space and storing
        the results."""
        # TODO Find smarter ways of extracting consensus.
        consensus = Seq(alignment.get_consensus())
        # If considering reverse primers, work on reverse complement of
        # consensus.
        if rev:
            consensus = consensus.reverse_complement()
            # Accounting for reversed positions.
            self._allowed_5p = reverse_inds(allowed_5p, len(consensus))
        else:
            self._allowed_5p = allowed_5p

        super().__init__(num_params_to_keep, consensus, alignment)

        self._consensus = str(consensus)
        self._allowed_lengths = allowed_lengths
        self._adapters = adapters
        self._rev = rev

    def find_best(self) -> None:
        """Finds the binding sequences with the lowest scores and stores them
        in <self._best_params>"""
        self._seq_iterator = self._get_iterator()

        # Binding and adapter always 5'-3'
        for binding in self._seq_iterator:
            param = self._seq_iterator.get_last_param()
            for i, adapter in enumerate(self._adapters):
                # Evaluate param score, if it belongs in the list, add it.
                param.set_adapter_num(i)
                self.score_param(param, binding, adapter)
                cond1 = param.get_score() < self._min_score
                cond2 = self._cur_num_params < self._num_params_to_keep
                if cond1 or cond2:
                    self._add_to_best_params(param)
                    # Get a new blank copy of the param.
                    param = self._seq_iterator.get_last_param()

    def get_best_params(self) -> List[BindingParams]:
        """Returns the best params found by this class."""
        if not self._rev:
            return self._best_params
        else:
            for param in self._best_params:
                param.rev_5p_ind(len(self._consensus))
            return self._best_params

    def _get_iterator(self) -> HomoSeqIterator:
        """Returns an iterator over all possible binding sequences."""
        return HomoSeqIterator(consensus=self._consensus,
                               allowed_5p=self._allowed_5p,
                               allowed_lengths=self._allowed_lengths)

    def score_param(self, param: BindingParams, binding_seq: str,
                    adapter_seq: str) -> None:
        """Returns the score of the given adapter-binding seq combo."""
        self._find_average_conservation(param, self._rev)
        dimer_score = get_dimer_score(binding_seq, adapter_seq)
        param.set_dimer_score(dimer_score)

    def run_get_best_params(self, num_children: int = 1) -> List[BindingParams]:
        """Iterates over all possible binding sequences and returns the best
        results. Iff num_children > 1, spawns <num_children> processes and
        splits the binding space between them."""
        if num_children > 1:
            return self._split_into_n(num_children)

    @classmethod
    def _split_class(cls, num: int):
        """Splits the class into <num> separate homo_seq scorers, each sampling
        the same number of binding sequences. To be used in threading."""
        raise NotImplementedError('This feature is not implemented.')
        # To avoid uneven distribution or responsibilities, divide into even
        # portions by splitting parameters evenly. Size and adapter number
        # may contribute to runtime, so avoid splitting those.
        # 1. 5' inds
        # 2. Adapter seqs
        # 3. Size


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

    _for_adapters:
        Allowable forward adapters, 5'-3'.

    _rev_adapters:
        Allowable reverse adapters, 5'-3'.

    """
    _consensus: str
    _rc_consensus: str

    _seq_iterator: HeteroSeqIterator

    _f_allowed_5p: List[int]
    _r_allowed_5p: List[int]

    _f_allowed_len: List[int]
    _r_allowed_len: List[int]

    _allowed_amp_len: List[int]
    _min_allowed_amp_len: int

    _for_adapters: List[str]
    _rev_adapters: List[str]

    def __init__(self, alignment: MSA, f_allowed_5p: List[int],
                 r_allowed_5p: List[int], f_allowed_len: List[int],
                 r_allowed_len: List[int], allowed_amp_len: List[int],
                 forward_adapters: List[str], reverse_adapters: List[str],
                 num_params_to_keep: int) -> None:
        """Initialises a HeteroSeqScorer that will evaluate all primer pairs
        over the given range."""
        consensus = Seq(alignment.get_consensus())
        super().__init__(num_params_to_keep, Seq(consensus),
                         alignment)

        consensus_cpy = Seq(str(consensus[:]))
        self._rc_consensus = str(consensus_cpy.reverse_complement())
        self._consensus = str(consensus_cpy)

        self._f_allowed_5p = f_allowed_5p
        self._r_allowed_5p = r_allowed_5p

        self._f_allowed_len = f_allowed_len
        self._r_allowed_len = r_allowed_len

        self._allowed_amp_len = allowed_amp_len
        self._min_allowed_amp_len = min(allowed_amp_len)

        self._for_adapters = forward_adapters
        self._rev_adapters = reverse_adapters

    def get_best_params(self) -> List[BindingPairParams]:
        """Returns the best params found by this class."""
        for param in self._best_params:
            param.r_params.rev_5p_ind(len(self._consensus))
        return self._best_params

    def find_best(self, V: bool = False) -> None:
        """Finds the binding sequences with the lowest scores and stores them
        in <self._best_params>"""
        self._seq_iterator = self._get_iterator()
        num = 0

        # Binding and adapter always 5'-3'
        for binding_pair in self._seq_iterator:
            fb, rb = binding_pair
            param = self._seq_iterator.get_last_param()
            for i in range(len(self._for_adapters)):
                num += 1
                if V:
                    if num % PRINT_EVERY == 0:
                        print(num, 'Possible pairs have been compared.')
                fa = self._for_adapters[i]
                ra = self._rev_adapters[i]
                # Compare binding-adapter pairs. Incorporate score if possible.
                self.score_param(param, fa, rb, ra, rb)
                cond1 = param.get_score() < self._min_score
                cond2 = self._cur_num_params < self._num_params_to_keep
                if cond1 or cond2:
                    param.f_params.set_adapter_num(i)
                    param.r_params.set_adapter_num(i)
                    self._add_to_best_params(param)
                    param = self._seq_iterator.get_last_param()

    def score_param(self, param: BindingPairParams, fa: str, fb: str, ra: str,
                    rb: str) -> None:
        """Returns the score of the given adapter-binding seq combo."""
        self._find_average_conservation(param.f_params, False, 'f')
        self._find_average_conservation(param.r_params, False, 'r')
        mean_cons = (param.f_params.get_mean_conservation() +
                     param.r_params.get_mean_conservation()) / 2
        dimer_score = get_heterodimer_score(fa, fb, ra, rb)
        param.set_mean_conservation(mean_cons)
        param.set_dimer_score(dimer_score)

    def _get_iterator(self) -> HeteroSeqIterator:
        """Returns the iterator that will be used to iterate over the allowable
        binding region combinations."""
        return HeteroSeqIterator(self._consensus, self._f_allowed_5p,
                                 self._r_allowed_5p, self._f_allowed_len,
                                 self._r_allowed_len, self._allowed_amp_len)

    @classmethod
    def _split_class(cls: C, num: int) -> List[C]:
        raise NotImplementedError('This feature is not implemented.')

    def run_get_best_params(self, num_children: int = 1) -> List[P]:
        pass


class BestPrimers:
    """A class intended to rigorously search for the best possible primers in
    a given region of an alignment.

    === Private Attributes ===

    _alignment:
        The alignment for which binding sequences will be analysed.

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

    _for_adapters:
        Allowable forward adapters, 5'-3'.

    _rev_adapters:
        Allowable reverse adapters, 5'-3'.

        _for_scorer: HomoSeqScorer
        Used to find the most performant
    _rev_scorer: HomoSeqScorer

    _pair_scorer: HeteroSeqScorer

    """
    _alignment: MSA

    _consensus: str
    _rc_consensus: str

    _seq_iterator: HeteroSeqIterator

    _f_allowed_5p: List[int]
    _r_allowed_5p: List[int]

    _f_allowed_len: List[int]
    _r_allowed_len: List[int]

    _allowed_amp_len: List[int]
    _min_allowed_amp_len: int

    _for_adapters: List[str]
    _rev_adapters: List[str]

    _for_scorer: HomoSeqScorer
    _rev_scorer: HomoSeqScorer
    _pair_scorer: HeteroSeqScorer

    _num_params_to_keep: int

    def __init__(self, alignment: MSA, f_allowed_5p: List[int],
                 r_allowed_5p: List[int], f_allowed_len: List[int],
                 r_allowed_len: List[int], allowed_amp_len: List[int],
                 forward_adapters: List[str], reverse_adapters: List[str],
                 num_params_to_keep: int) -> None:
        """Initialises a HeteroSeqScorer that will evaluate all primer pairs
        over the given range."""

        self._alignment = alignment

        self._num_params_to_keep = num_params_to_keep
        consensus = Seq(alignment.get_consensus())

        consensus_cpy = Seq(str(consensus[:]))
        self._rc_consensus = str(consensus_cpy.reverse_complement())
        self._consensus = str(consensus_cpy)

        self._f_allowed_5p = f_allowed_5p
        self._r_allowed_5p = r_allowed_5p

        self._f_allowed_len = f_allowed_len
        self._r_allowed_len = r_allowed_len

        self._allowed_amp_len = allowed_amp_len
        self._min_allowed_amp_len = min(allowed_amp_len)

        self._for_adapters = forward_adapters
        self._rev_adapters = reverse_adapters

        self._init_scorers()


    def _init_scorers(self) -> None:
        """Initialises the scorers."""
        self._for_scorer = HomoSeqScorer(self._alignment, self._f_allowed_5p,
                                         self._f_allowed_len,
                                         self._for_adapters,
                                         self._num_params_to_keep, False)
        self._rev_scorer = HomoSeqScorer(self._alignment, self._r_allowed_5p,
                                         self._r_allowed_len,
                                         self._rev_adapters,
                                         self._num_params_to_keep, True)
        self._pair_scorer = HeteroSeqScorer(self._alignment, self._f_allowed_5p,
                                            self._r_allowed_5p,
                                            self._f_allowed_len,
                                            self._r_allowed_len,
                                            self._allowed_amp_len,
                                            self._for_adapters,
                                            self._rev_adapters,
                                            self._num_params_to_keep)

    def get_n_best(self, n: int, V: bool = False) -> List[BindingPairParams]:
        """Returns a list of lowest scoring binding pair params."""
        self._for_scorer.find_best()
        self._rev_scorer.find_best()
        if V:
            print('Completed individual analysis, moving onto pairwise.')
        self._pair_scorer.find_best(V)

        def find_in(bind_param: BindingParams, params: List[BindingParams]) \
                -> Union[bool, BindingParams]:
            """Returns the bind_param if it is present in params, else returns
            false."""
            for param in params:
                if param == bind_param:
                    return param
            return False

        min_for = self._for_scorer.get_best_params()
        min_rev = self._rev_scorer.get_best_params()
        min_pair = self._pair_scorer.get_best_params()

        pair_params_to_bind_params = {}

        # If the individual params for any pair of params have been found,
        # store them.
        for pair_param in min_pair:
            found_params = [find_in(pair_param.f_params, min_for),
                            find_in(pair_param.r_params, min_rev)]
            if not found_params[1]:
                found_params.pop(1)
            elif not found_params[0]:
                found_params.pop(0)
            pair_params_to_bind_params[pair_param] = found_params

        # Dictionaries for pairs where both, one and none of the pair were
        # found in the homo minimum lists.
        both = {}
        one = {}
        none = {}

        # Calculate score for each dimer in each case.
        for pair_param in min_pair:
            if len(pair_params_to_bind_params[pair_param]) == 2:
                f_param = pair_params_to_bind_params[pair_param][0]
                r_param = pair_params_to_bind_params[pair_param][1]
                both[pair_param] = f_param.get_dimer_score() ** DIMER_WEIGHT * \
                                   r_param.get_dimer_score() ** DIMER_WEIGHT
                both[pair_param] = both[pair_param] * pair_param.get_score()
                pair_param.set_final_score(both[pair_param])


            elif len(pair_params_to_bind_params[pair_param]) == 1:
                param = pair_params_to_bind_params[pair_param][0]
                one[pair_param] = param.get_dimer_score() ** DIMER_WEIGHT
                one[pair_param] = one[pair_param] * pair_param.get_score()
                pair_param.set_final_score(both[pair_param])
            else:
                none[pair_param] = pair_param.get_score()
                pair_param.set_final_score(both[pair_param])

        best_both = list(both.keys())
        find_score = lambda a: both[a]
        best_both.sort(key=find_score)

        best_one = list(one.keys())
        find_score = lambda a: one[a]
        best_one.sort(key=find_score)

        best_none = list(none.keys())
        find_score = lambda a: none[a]
        best_none.sort(key=find_score)

        best_params = best_both + best_one + best_none

        to_rtrn = best_params[:min(n, len(best_params))]


        for param in to_rtrn:
            fb = binding_param_to_seq(self._consensus, param.f_params,
                                      rev=False)
            fa = self._for_adapters[param.f_params.get_adapter_num()]
            self._for_scorer.score_param(param.f_params, fb, fa)
            rb = binding_param_to_seq(self._consensus, param.r_params,
                                      rev=True)
            ra = self._rev_adapters[param.r_params.get_adapter_num()]
            self._rev_scorer.score_param(param.r_params, rb, ra)

        return best_params[:min(n, len(best_params))]


def vis_score(bps: BindingPairParams, ind: str = '    ') -> str:
    """Returns a string containing the attributes of the given binding
    parameters. Indents all lines with indent."""
    fp = bps.f_params
    rp = bps.r_params
    wcs = bps.get_mean_conservation() ** CONS_WEIGHT
    wds = (bps.get_dimer_score() ** DIMER_WEIGHT *
                            fp.get_dimer_score() ** DIMER_WEIGHT *
                            rp.get_dimer_score() ** DIMER_WEIGHT)
    return ''.join([
        ind, 'Final Score: ', str(bps.get_final_score()), '\n',
        ind, 'Weighted Conservation Score: ', str(wcs), '\n',
        ind, 'Weighted Dimer Score: ', str(wds), '\n',
        ind, 'Pair Dimer Score: ', str(bps.get_dimer_score()), '\n',
        ind, 'Forward Dimer Score: ', str(fp.get_dimer_score()), '\n',
        ind, 'Reverse Dimer Score: ', str(rp.get_dimer_score()), '\n',
        ind, 'Average Conservation: ', str(bps.get_mean_conservation()), '\n',
        ind, 'Forward Binding Range: ', str(fp.get_5p_ind()), '-',
        str(fp.get_5p_ind() + fp.get_len() - 1), '\n',
        ind, 'Reverse Binding Range: ', str(rp.get_5p_ind() - rp.get_len() + 1), '-',
        str(rp.get_5p_ind()), '\n',
        ind, 'Forward Binding Size: ', str(fp.get_len()), '\n',
        ind, 'Reverse Binding Size: ', str(rp.get_len()), '\n',
        ind, 'Forward Adapter Number: ', str(fp.get_adapter_num()), '\n',
        ind, 'Reverse Adapter Number: ', str(fp.get_adapter_num()), '\n',
        ])

def get_seqs(bps: BindingPairParams, consensus: str, ind: str = '    ') -> str:
    """Returns a text representaition of the sequences specified by bps."""
    fp = bps.f_params
    rp = bps.r_params
    return ''.join([
        ind, 'Forward Binding Range:', str(fp.get_5p_ind()), '-',
        str(fp.get_5p_ind() + fp.get_len() - 1),
        '\n Forward Binding Sequence \n',
        ind, '5\' ', binding_param_to_seq(consensus, fp, rev=False), ' 3\'\n',
        ind, 'Reverse Binding Range:', str(rp.get_5p_ind() - rp.get_len() + 1),
        '-',
        str(rp.get_5p_ind()), '\n',
        'Reverse Binding Sequence \n',
        ind, '5\' ',
        str(Seq(binding_param_to_seq(consensus, rp, rev=True)).complement()),
        ' 3\'\n',
    ])
