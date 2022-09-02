from test_files.multiplex_testing.test_generate_spacer_seqs import rand_seq
import random

def test_basic() -> None:
    t1 = primer3_adapter_float.calc_heterodimer_score(
        ('G' + 'A' * 90 + 'T').upper(), ('C' + 'A' * 90 + 'TT').upper())
    t2 = primer3_adapter_float.calc_heterodimer_score(
        'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAT',
        'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAATT')
    print(t1, t2)


def test_basic_max_len() -> None:
    t1 = primer3_adapter_float.calc_heterodimer_score(
        'GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTCATCTTGCAGTGGCGGAACTGCTTGT',
        'ACACTCTTTCCCTACACGACGCTCTTCCGATCTTGTTTTCTCCTACTTGGCTGCAATCTGGAAGG')
    print(t1)


def test_basic_max_len_iter() -> None:
    min_seq_len = 1
    max_seq_len = 200
    for _ in range(1000):
        seq = rand_seq(random.randint(min_seq_len, max_seq_len))
        seq_len = len(seq)
        trimmed_seq = seq_len_check(seq)
        if seq_len <= 60:
            assert seq == trimmed_seq
        else:
            assert seq[-60:] == trimmed_seq

def test_basic_max_len_iter_double() -> None:
    min_seq_len = 30
    max_seq_len = 200
    for _ in range(1000):
        seq1 = rand_seq(random.randint(min_seq_len, max_seq_len))
        seq2 = rand_seq(random.randint(min_seq_len, max_seq_len))
        trimmed_seq1, trimmed_seq2 = seq_len_check_2(seq1, seq2)
        if len(seq1) <= 60:
            assert seq1 == trimmed_seq1
        else:
            assert seq1[-60:] == trimmed_seq1
            seq1 = seq1[-60:]

        if len(seq2) <= 60:
            assert seq2 == trimmed_seq2
        else:
            assert seq2[-60:] == trimmed_seq2
            seq2 = seq2[-60:]

        t1 = primer3_adapter_float.calc_heterodimer_score(seq1, seq2)

        t2 = primer3_adapter_float.calc_heterodimer_score(trimmed_seq1,
                                                          trimmed_seq2)
        assert t1 == t2


