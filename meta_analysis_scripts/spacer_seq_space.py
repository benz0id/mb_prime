import threading
from statistics import stdev
from hetero_spacer_generator import hetero_spacer_generator as hsg
from test_files.fixtures_and_helpers import gen_random_seq
from math import factorial
from typing import Tuple, List, Dict

'''
A script to estimate the size sequence space of the average heterogeneity
spacer alignment.
'''
# Verbosity
V = True
# Number of sequences to sample
N = 1000000
# Some number that divides 100
alert_every_prcnt = 5
# Parameter for the number of spacers to be generated.
num_hetero_bases = 12
max_spacer_length = 12
seq_len = 12
GC = 0.41


def main():
    run()


def nPr(n: int, r: int) -> int:
    '''Returns <n> permute <r>. Precondition n >= r.
    '''
    return int(factorial(n) / factorial(n - r))


def align_to_space(alignment: Tuple[int, int, int, int]) -> int:
    """Given some <alignment>, returns the number of possible sequences that
    could occupy that alignment.
    Precondition: <alignment> must be sorted LTG."""
    num_permutations = [nPr(4, 4) ** (alignment[0]),
                        (nPr(3, 3) ** (alignment[1] - alignment[0])),
                        (nPr(2, 2) ** (alignment[2] - alignment[1]))]
    # Calculate product of all of these permutations.
    prod = 1
    for i in range(len(num_permutations)):
        prod = prod * max(num_permutations[i], 1)
    return prod


def sample_n_spacers(n: int,
                     spacers_dict: Dict[int, List[float]] or
                                   List[float]) -> None:
    # Setting up required tools.
    seq_spaces = []
    global do_threading
    SAG = hsg.SpacerAlignmentGen(max_spacer_length, num_hetero_bases)
    name = threading.current_thread().name
    # Generate spacers for <N> randoms sequences.
    for i in range(int(100 / alert_every_prcnt)):
        for j in range(int(n / (int(100 / alert_every_prcnt)))):
            seq = gen_random_seq(12, GC)
            spacer_aligns = SAG.get_all_spacer_combos(seq)
            SAG.sort_spacer_combos(spacer_aligns)
            best_spacer_align = spacer_aligns[0]
            ats = align_to_space(best_spacer_align)
            seq_spaces.append(ats)
        if V:
            print(name + ' is ' + str(
                (i + 1) * alert_every_prcnt) + "% complete.")
    print(len(seq_spaces))
    if isinstance(spacers_dict, dict):
        spacers_dict[sum(spacers_dict.keys())] = seq_spaces
    else:
        spacers_dict.extend(seq_spaces)
    if V:
        print(name + " complete")


def run() -> int:
    seq_spacers = []
    sample_n_spacers(N, seq_spacers)
    if V:
        print(
            "NUM NON-THREADING SPACERS CHECKED = " + str(len(seq_spacers)))
    avg = sum(seq_spacers) / len(seq_spacers)
    print("The computed average sequence space is " +
          str(avg) +
          ", with a standard deviation of " +
          str(int(stdev(seq_spacers))))
    return avg



def calc_for_gc():
    global GC
    n_steps = 100
    GC_to_ss = {}
    for GC20 in range(1, n_steps):
        GC = GC20 / n_steps
        if V:
            print("Calculating for GC = " + str(GC))
        GC_to_ss[GC] = run()
    print(GC_to_ss)


if __name__ == '__main__':
    main()
