from collections import Iterable, Iterator
from abc import ABC, abstractmethod
from typing import List
from Bio.Seq import Seq
from multiprocessing import Process
from bisect import bisect_right


# Where to change if new parameters are found:
#   Binding Params
#
#
#


class Param:
    """Stores some parameters and their associated overall score.

    _score:
        The score of this parameter"""

    _score: float

    def __init__(self) -> None:
        self._score = 0

    def __lt__(self, other) -> bool:
        return self._score < other.get_score()

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

    def __init__(self, fivep_ind: int, length: int, adapter_num: int) -> None:
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

    def get_adapter_num(self) -> None:
        """Returns this param's adapter number."""
        print(self._adapter_num)


class BindingPairParams(Param):
    """The parameters used to construct a pair of forward and reverse binding
    sequences.

    === Private Attributes ===
    _f_params:
        Parameters of the forward binding seq.
    _r_params:
        Parameters of the reverse binding seq."""

    _f_params: BindingParams
    _r_params: BindingParams

    def __init__(self, f_params: BindingParams, r_params: BindingParams) \
            -> None:
        """Initialises this class with the given foward and reverse binding
        parameters."""
        super().__init__()
        self._f_params = f_params
        self._r_params = r_params


class BindingSeqIterator(Iterable, Iterator, ABC):
    """A class that iterates over all possible permutations of a given set of
    parameters of binding sequences, returning the binding sequences specified
    by those parameters.
    """

    @abstractmethod
    def get_param(self) -> Param:
        """Returns the parameters of the binding sequence returned by the last
        call as a Param object"""
        pass


class BindingSeqScorer(ABC):
    """A class designed to score and collect the best outputs of a
    BindingSeqIterator.

    === Private Attributes ===
    _best_params:
        A list containing the highest scoring parameters. Sorted G-L
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

    _best_params: List[Param]
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

    def _add_to_best_params(self, param: Param) -> None:
        """Adds <param> to <self._best_params>, maintaining G-L order."""
        ind = bisect_right(self._best_params, param)
        self._best_params.insert(ind, param)
        # Remove last ind if desired length has been reached.
        if len(self._best_params) > self._num_params_to_keep:
            self._best_params.pop(-1)
        else:
            self._cur_num_params += 1

    def _add_if_better(self, score: float) -> None:
        """Adds the parameters used to construct the last binding sequence(s) to
        self._best_params iff it is better than the worst score in the list or
        the list is incomplete."""
        if score > self._min_score or \
                self._cur_num_params < self._num_params_to_keep:
            param = self._seq_iterator.get_param()
            self._add_to_best_params(param)

    def _split_into_n_(self, num_children: int):
        """Finds and returns the parameters that produce the best binding
        sequences by splitting the task among <num_threads> children."""
        split_classes = self._split_class(num_children)
        children = []
        for bss in split_classes:
            children.append(Process(target=bss._find_best))

    @abstractmethod
    def get_best_params(self, num_children: int = 1) -> List[Param]:
        """Return the parameters that produce the best binding sequences."""
        pass

    @abstractmethod
    def _get_iterator(self) -> BindingSeqIterator:
        """Return the iterator from which this BindingSeqScorer will extract
        binding sequences."""
        pass

    @abstractmethod
    def _split_class(self, num: int) -> List[BindingSeqScorer]:
        """Splits the binding sequence space into <num> other BindingSeqScorers
        and returns them."""
        pass

    @abstractmethod
    def _find_best(self) -> None:
        """Completely iterates through <self._seq_iterator> and stores the best
        scoring params in <self._best_params>"""
        pass


class HomoSeqScorer:
    """Evaluates the features of primers produced by all possible
    permutations of some construction parameters.

    for_binding |       5'-GCATGGTGATCGT-3'
    <consensus> | 5'-ATCGATGCATGGTGATCGTAGGTCGGTAGTCGTCAGTCGTGTAGTCGG-3'
    rev_binding |                              3'-AGCAGTCAGCACA-5'

    === Private Attributes ===

    _alignment:
        The alignment for which binding sequences will be analysed.

    _rev:
        Iff true will consider primers for the reverse compliment of
        alignment consensus sequence.

    _max_len:
        The maximum length of any binding region to be considered.
    _min_len:
        The minimum length of any binding region to be considered.

    _region_start:
        The index of the last base in region of adapter to which this primer
        can bind.
    _region_end:
        The index of the last base in region of adapter to which this primer
        can bind.

    _adapters:
        A list of adapter sequences 5' - 3'.

    _scores_array: numpy array. Stores the scores of all possible primers.

        _scores_array[primer_5p_ind][primer_len][adapter_number] == score
            for that primer.

    """


class HeteroSeqScorer:
    """
    Evaluates the features of each allowable pair of forward and reverse
    binding sequence.

    _alignment:

    _f_region_start:
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
