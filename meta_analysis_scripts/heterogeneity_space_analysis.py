from pathlib import Path
from typing import List, Tuple, Union
import matplotlib.pyplot as plt
from Bio.Seq import Seq
import winsound
from hetero_spacer_generator.get_random_seqs import gen_hetero_set
from hetero_spacer_generator.primer_tools import MBPrimerBuilder
from hetero_spacer_generator.hetero_spacer_generator import \
    HeteroGen
from numpy import std

PICKLE_PATH = Path('C:\\Users\\bfern\\PycharmProjects\\mb_prime\\pickle')
FONT_SIZE = 18
SPACER_LENGTH = 12
NUM_HOMO = 10000
NUM_HETERO = 40
BT = 500
COLOUR = 'r'
OUT_PATH = Path(
    'C:\\Users\\bfern\\OneDrive - University of Toronto\\Desktop\\Good Sets')
SET_NAME = 'Good'

fa1 = Seq('ACACTCTTTCCCTACACGACGCTCTTCCGATCT')
ra1 = Seq('GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT')
fb1 = Seq("GTCATCTTCTTCTGCTACG")
rb1 = Seq("GCCGGAATGGTCATGAAGA")

regular_set = fa1, ra1, fb1, rb1

fa2 = Seq('ACACTCTTTCCCTACACGACGCTCTTCCGATCT')
ra2 = Seq('GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT')
fb2 = Seq("GTCATCTTCTTCTGAGATC")
rb2 = Seq("GCCGGAATGGTCATAGATC")

bad_set = fa2, ra2, fb2, rb2

plt.rcParams.update({'font.size': FONT_SIZE})

def main():
    # HuRho_1 binding and adapter sequences.
    fa, ra, fb, rb = regular_set
    winsound.Beep(1000, BT)
    # Construct incomplete primers
    fp = MBPrimerBuilder(adapter_seq=fa, binding_seq=fb)
    rp = MBPrimerBuilder(adapter_seq=ra, binding_seq=rb)

    # Get and display scores.
    fhs, rhs = get_homogeneity_scores(NUM_HOMO, fp, rp)
    winsound.Beep(1000, int(BT / 2))
    winsound.Beep(1200, int(BT / 2))
    hetero_scores = get_heterogeneity_scores(NUM_HETERO, fp, rp)
    mean = lambda a: sum(a) / len(a)
    get_scores_histogram(fhs, title='Forward Score Space')
    get_scores_histogram(rhs, title='Reverse Score Space')
    get_scores_histogram(hetero_scores, title='Pair Score Space')
    print(mean(fhs))
    print(mean(rhs))
    print(mean(hetero_scores))
    print(std(fhs))
    print(std(rhs))
    print(std(hetero_scores))
    winsound.Beep(1000, int(BT / 3))
    winsound.Beep(1200, int(BT / 3))
    winsound.Beep(1000, int(BT / 3))

    """GOOD:
    37.9849
38.6374
33.754375

BAD
40.6406
40.3225
33.710625"""


# Script to analyse the score space of metabarcoding primers.
def get_homogeneity_scores(n: int, incomplete_forward: MBPrimerBuilder,
                           incomplete_reverse: MBPrimerBuilder) \
        -> Tuple[List[int], List[int]]:
    """Generates n random spacer sets, for the given incomplete forward and
    reverse primers, returns the scores as calculated by a pairwise sorter.
    Returns scores forward first, and reverse second."""

    # Generate and extract scoring objects.
    hg = HeteroGen(rigour=100,
                   num_hetero=SPACER_LENGTH,
                   max_spacer_length=SPACER_LENGTH)
    rg = hg.get_primer_gen()
    sf = rg.get_spacer_sorter()

    # Select best spacer combos.
    f_spacer = hg.get_all_spacer_combos(
        incomplete_forward.get_binding_seq())[0]

    r_spacer = hg.get_all_spacer_combos(
        incomplete_reverse.get_binding_seq())[0]

    # Generate Random Spacers.
    forward_spacer_seqs = gen_hetero_set(incomplete_forward,
                                         f_spacer,
                                         n)
    reverse_spacer_seqs = gen_hetero_set(incomplete_reverse,
                                         r_spacer,
                                         n)

    # Compute spacer scores.
    sf._build_full(SPACER_LENGTH, SPACER_LENGTH, incomplete_forward,
                   incomplete_reverse, forward_spacer_seqs, reverse_spacer_seqs)

    sf._evaluate_scores_single()

    # Extract scores
    f_scores = sf._for_single_scores
    r_scores = sf._rev_single_scores

    return f_scores, r_scores


def get_heterogeneity_scores(n: int, incomplete_forward: MBPrimerBuilder,
                             incomplete_reverse: MBPrimerBuilder) -> List[
    float]:
    """Generates n random spacer sets, for the given incomplete forward and
    reverse primers, returns the scores as calculated by a pairwise sorter.
    Returns scores forward first, and reverse second."""

    # Generate and extract scoring objects.
    hg = HeteroGen(rigour=100,
                   num_hetero=SPACER_LENGTH,
                   max_spacer_length=SPACER_LENGTH)
    rg = hg.get_primer_gen()
    sf = rg.get_spacer_sorter()

    # Select best spacer combos.
    f_spacer = hg.get_all_spacer_combos(
        incomplete_forward.get_binding_seq())[0]

    r_spacer = hg.get_all_spacer_combos(
        incomplete_reverse.get_binding_seq())[0]

    # Generate Random Spacers.
    forward_spacer_seqs = gen_hetero_set(incomplete_forward,
                                         f_spacer,
                                         n)
    reverse_spacer_seqs = gen_hetero_set(incomplete_reverse,
                                         r_spacer,
                                         n)

    # Generate matrix of spacer pair combinations.
    sf._build_full(SPACER_LENGTH, SPACER_LENGTH, incomplete_forward,
                   incomplete_reverse, forward_spacer_seqs, reverse_spacer_seqs)
    sf._num_pairings_to_compare = n
    sf._construct_half_sets()
    sf._construct_full_set_matrix()
    sets = sf._fullsets

    # Score each set.
    scores = []
    for row in sets:
        for set in row:
            set.apply_criteria(sf._pair_criteria,
                               sf._pair_criteria_weights_list)
            scores.append(set.get_score())

    return scores


def get_scores_histogram(x: List[Union[int, float]],
                         title: str = 'Distribution of Set Scores'
                         , c: str = COLOUR) -> None:
    """Generates a scores histogram given some scores in x."""
    num_bins = len(set(x))
    plt.figure(figsize=(9, 7), dpi=80)
    plt.hist(x, num_bins, facecolor=c)
    plt.xlabel('Score')
    plt.ylabel('Number of Sets')
    plt.title(title)
    plt.grid(True)
    plt.savefig(OUT_PATH / (title + ' ' + SET_NAME))
    plt.close()


if __name__ == '__main__':
    main()
