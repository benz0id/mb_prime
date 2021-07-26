from typing import Tuple, List, Callable
from primer_tools import MBPrimerBuilder, MBPrimer
from Bio.Seq import Seq


# A set of four alignments, i.e. the amount by which each primer should be
# shifted to the right in a set of four to ensure heterogeneity.
SpacerAlignment = Tuple[int, int, int, int]

# A list of four sets of forward or reverse heterogeneity spacer sequences.
SpacerSet = Tuple[Seq, Seq, Seq, Seq]

# Where <SpacersSets> is the set of all spacers seqs, <MBPrimerBuilder> is the
# incomplete spacer. Returns a list of all scores.
PairWiseCriterionSingle = Callable[[SpacerSet, MBPrimerBuilder], List[int]]
# Similar to PairWiseCriterionSingle except that it evaluates traits between the
# two sets and primers and returns a matrix. Forward primers come first both in
# indexing the matrix and calling the functions.
SimpleCriterion = Callable[[Seq, Seq], int]
