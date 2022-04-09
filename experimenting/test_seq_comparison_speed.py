from typing import List, Tuple
from Bio.Seq import Seq
import random
import timeit
from hetero_spacer_generator.sequence_tools import SeqAnalyzer
from statistics import stdev

# What is use throughout the program to analyse sequences.
sa = SeqAnalyzer(degen=False)
random.seed(a=123456789, version=1)

# Length of the two sequences to be compared
seq1_len = 33
seq2_len = 33

# The number of different randomly generated sequences to compare.
num_rand_seqs = 100
# The number of times to re-run each timing. Increases accuracy.
num_repetitions = 100
# A factor by which to multiple the final time. To be used when expecting to
# do more than one string comparison.
factor = 1


def rand_seq(num: int) -> str:
    """Returns a random nucleotide sequence with <num> bases."""
    bases = ['A', 'T', 'C', 'G']
    seq_lst = [random.choice(bases) for _ in range(num)]
    return ''.join(seq_lst)


def base_comparison_times(seq1: Seq, seq2: Seq,
                          num_reps: int = num_repetitions) \
        -> Tuple[float, float]:
    """Returns the average times to find the consecutive and total
    complementarity between <seq1> 5' - 3' and <seq2> 3' - 5'. Performs num_reps
    repetitions."""
    get_consec = sa.get_consec_complementarity
    get_total = sa.get_non_consec_complementarity
    find_consec = lambda: sa.comp_seqs_any_overlap(seq1, seq2[::-1], get_consec)
    find_total = lambda: sa.comp_seqs_any_overlap(seq1, seq2[::-1], get_total)
    consec_time = timeit.timeit(find_consec, number=num_reps) / num_reps
    total_time = timeit.timeit(find_total, number=num_reps) / num_reps
    return consec_time, total_time


def time_for_n_random(n: int = num_rand_seqs) -> Tuple[float, float, float]:
    """Returns the average time to find the consecutive and total
    complementarity between two randomly generated sequences as well as the
    standard deviation between runs. Will repeat <n> times."""
    consec_times = []
    total_times = []
    for _ in range(n):
        seq1 = rand_seq(seq1_len)
        seq2 = rand_seq(seq2_len)
        ct, tt = base_comparison_times(Seq(seq1), Seq(seq2))
        consec_times.append(ct)
        total_times.append(tt)
    return sum(consec_times) / n, sum(total_times) / n, \
           stdev(consec_times) + stdev(total_times)


def main():
    consec_time, total_time, stdv = time_for_n_random()
    print("Total Complementarity", total_time)
    print("Continuous Complementarity", consec_time)
    print("Combined, Factored Time:", (total_time + consec_time) * factor,
          ' +/-', stdv * factor)


if __name__ == "__main__":
    main()
