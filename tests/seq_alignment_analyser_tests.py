from seq_alignment_analyser.seq_alignment_analyser import SeqAlignAnalyser
from Bio.Seq import Seq
from fixtures_and_helpers import gen_random_seq, gen_random_seq_degen, \
    SMALL_SAMPLE_SIZE, gen_seq_list
import random



def test_degen_check() -> None:

    saa = SeqAlignAnalyser()
    for i in range(SMALL_SAMPLE_SIZE):
        saa.set_seqs(gen_seq_list(5, 50, variance=5, degen=False))
        saa._do_degen_check()
        assert saa._base_matches_generic == saa._base_matching_method

    for i in range(SMALL_SAMPLE_SIZE):
        saa.set_seqs(gen_seq_list(5, 50, variance=5, degen=True))
        saa._do_degen_check()
        assert saa._base_matches_generic != saa._base_matching_method
        for i in range(len(saa._seqs)):
            assert len(saa._seqs_degen_array) == len(saa._seqs)


def test_gen_complements() -> None:
    saa = SeqAlignAnalyser()
    for i in range(SMALL_SAMPLE_SIZE):
        num_seqs = random.randrange(0, SMALL_SAMPLE_SIZE)
        saa.set_seqs(gen_seq_list(5, 50, variance=5, degen=False))
        saa._gen_complements()
        assert len(saa._seqs_complement) == len(saa._seqs)

    for i in range(SMALL_SAMPLE_SIZE):
        num_seqs = random.randrange(0, SMALL_SAMPLE_SIZE)
        saa.set_seqs(gen_seq_list(5, 50, variance=5, degen=True))
        saa._gen_complements()
        assert len(saa._seqs_complement) == len(saa._seqs)


seqs_test_hga_1 = [Seq('AAAA'),
                   Seq('AAAT'),
                   Seq('AATC'),
                   Seq('ATTG')]


def test_gen_heterogeneity_array_same_len() -> None:
    saa = SeqAlignAnalyser()
    saa.set_seqs(seqs_test_hga_1)
    saa._do_degen_check()
    saa._gen_complements()
    saa._gen_heterogeneity_array()
    for i in range(0, 3):
        for j in range(i + 1, 3):
            assert saa._bases_match(0, i, j)

    for i in range(0, 3):
        for j in range(i + 1, 3):
            assert not saa._bases_match(3, i, j)

seqs_test_hga_2 = [Seq('MAAA'),
                   Seq('AATTAA'),
                   Seq('ATTCGSM'),
                   Seq('NTTGNNN')]


def test_gen_heterogeneity_array_dif_len_degen() -> None:
    saa = SeqAlignAnalyser()
    saa.set_seqs(seqs_test_hga_2)
    saa._do_degen_check()
    saa._gen_complements()
    saa._gen_heterogeneity_array()
    for i in range(0, 3):
        for j in range(i + 1, 3):
            assert saa._bases_match(0, i, j)

    for i in range(0, 3):
        for j in range(i + 1, 3):
            assert not saa._bases_match(3, i, j)


seqs_test_nm_1 = [Seq('AAAAAA'),
                  Seq('AAATTT'),
                  Seq('AATTCC'),
                  Seq('ATTTCG')]


def test_num_matching() -> None:
    saa = SeqAlignAnalyser()
    saa.set_seqs(seqs_test_nm_1)
    saa._gen_heterogeneity_array()
    saa._gen_num_matching()

    assert saa._seqs_num_matching[0] == 6
    assert saa._seqs_num_matching[1] == 3
    assert saa._seqs_num_matching[2] == 2
    assert saa._seqs_num_matching[3] == 3
    assert saa._seqs_num_matching[4] == 1
    assert saa._seqs_num_matching[5] == 0







