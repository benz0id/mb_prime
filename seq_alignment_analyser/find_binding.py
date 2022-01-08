# Classes, functions to find ideal binding sequences.
from Bio.Seq import Seq
from typing import List, Tuple

# A highly rigorous search for the best binding sequence in a region with the
# given parameters.


def find_best_binding_pair(region: Seq, for_adapters: List[Seq],
                           rev_adapters: List[Seq], for_size: Tuple[int, int],
                           rev_size: Tuple[int, int], min_amp_size: int,
                           max_amp_size: int) -> Tuple[Seq, Seq]:
    """Finds the best forward and reverse binding sequences of lengths
    <for_size> and <rev_size> respectfully.

    Returns tuple (for_binding, rev_binding).

    for_binding |       5'-GCATGGTGATCGT-3'
    <region>    | 5'-ATCGATGCATGGTGATCGTAGGTCGGTAGTCGTCAGTCGTGTAGTCGG-3'
    rev_binding |                              3'-AGCAGTCAGCACA-5'

    Directionality:
        region: 5' - 3' sense
        for_adapters: 5' - 3' sense
        rev_adapters: 5' - 3' antisense
    """

    # for_binding |       5'-GCATGGTGATCGT-3'
    # for_region  | 5'-ATCGATGCATGGTGATCGTAGGTCGGTAGTCGTCAGTCGTGTAGTCGG-3'
    for_region = region
    # rev_binding |       5'-ACACGACTGACGA-3'
    # rev_region  | 5'-CCGACTACACGACTGACGACTACCGACCTACGATCACCATGCATCGAT-3'
    rev_region = region.reverse_complement()

    # ind: the first base in each binding sequence with respect to it's region
    # (i.e. the 5' end).
    max_ind = len(region) - min_amp_size

    # for_scores[primer_ind][adapter_num] = score of that binding seq
    for_scores = gen_binding_scores(region, for_adapters, for_size, max_ind)
    rev_scores = gen_binding_scores(rev_region, rev_adapters, rev_size, max_ind)


def gen_binding_scores(region: Seq, adapters: List[Seq], binding_size: int,
                       max_ind: int) -> List[List[int]]:
    """Generates an array of scores for each binding sequence of length
    <binding_size> in region. Only indices of region less than or equal to
     <max_ind> will be considered as valid binding sequence sites. Scores are
     indicative of each primers inherent homodimerisation potential.

    # binding_seq |       5'-GCATGGTGATCGT-3'
    # region      | 5'-ATCGATGCATGGTGATCGTAGGTCGGTAGTCGTCAGTCGTGTAGTCGG-3'

    return[ind][adapter_num] == score
    """

    scores = []
    # For each allowable index, create a potential binding sequence.
    for i in range(max_ind + 1):
        scores.append([])
        binding_seq = region[i:i + binding_size]

        for adapter_tup in enumerate(adapters):



