from typing import List, Tuple

from hetero_spacer_generator.primer_tools import PairwisePrimerSet, PrimerSet
from meta_tools.analysis_tools import Heterodimer, Homodimer, Parameter, Result




def get_most_stable_pairwise(primer_set: PrimerSet) -> \
    Tuple[List[Homodimer], List[Heterodimer]]:
    """Returns the most stable Homo and Heterodimers that might be formed
    when using this set."""
