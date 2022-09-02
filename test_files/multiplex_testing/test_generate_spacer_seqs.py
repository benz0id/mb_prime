from src.multiplex_spacer_generator.find_best_spacer_combo import FindSpacerCombo
from multiplex_spacer_generator.generate_spacer_seqs import *
import random
from multiprocessing import Queue
from test_files.fixtures_and_helpers import add_stdout_handler, \
    standard_pytest_run, configure_log_out

NUM_THREADS = 8

add_stdout_handler()

configure_log_out('hetero_seq_generation')


def test_hetero_seq_generation_threaded() -> None:
    max_het = 12
    max_seqs = 12
    num_reps = 20
    max_adapter_len = 24
    max_seq_len = 24
    spacer_finder_runtime = 10
    seq_generation_runtime = 10

    for _ in range(num_reps):
        het_len = random.randint(6, max_het)
        num_seqs = random.randint(6, max_seqs)
        f_binding = [rand_seq(max_seq_len) for _ in range(num_seqs)]
        r_binding = [rand_seq(max_seq_len) for _ in range(num_seqs)]
        f_adapters = [rand_seq(max_adapter_len) for _ in range(num_seqs)]
        r_adapters = [rand_seq(max_adapter_len) for _ in range(num_seqs)]
        num_structs_to_view = num_seqs ** 2 // 2

        spacer_finder = FindSpacerCombo(spacer_finder_runtime, NUM_THREADS,
                                        f_binding, het_len)
        f_spacers = spacer_finder.run()
        spacer_finder = FindSpacerCombo(spacer_finder_runtime, NUM_THREADS,
                                        r_binding, het_len)
        r_spacers = spacer_finder.run()

        final_pool = get_best_heterogeneity_spacer_seqs_threadable(
            f_adapters, f_binding, r_adapters, r_binding, f_spacers, r_spacers,
            seq_generation_runtime, num_structs_to_view, NUM_THREADS)
        assert final_pool






def test_gen_heterogeneity_alignment() -> None:
    seqs = [
        'ATGCTAGCTGACTGATCGAGCATGACT',
        'CAGCTAGCTAGTCGATCGACTGAGCAT',
        'ACTAGTCGATGCTGCATGCTAGTAGCT',
        'ATGCTAGCTGACTGATCGATGGTAGCT',
        'TAGCTGTAGCTAGCTACTGACTGATCG'
    ]

    het_spacers = [0, 1, 2, 3, 4]

    arr = gen_heterogeneity_alignment(seqs, het_spacers)


def test_get_rand_hetero_seqs() -> None:
    seqs = [
        'ATGCTA',
        'CAGCTA',
        'ACTAGT',
        'ATGCTA',
        'TAGCTG'
    ]

    het_spacers = [0, 1, 2, 3, 4]

    spacer_seqs = get_rand_hetero_seqs(seqs, het_spacers)

    for i, spacer in enumerate(spacer_seqs):
        assert len(spacer) == het_spacers[i]


rand_seq = lambda a: ''.join([random.choice('ATCG') for _ in range(a)])


def test_get_rand_hetero_seqs_random_seqs() -> None:
    max_het = 20
    max_seqs = 20
    num_reps = 10

    for _ in range(num_reps):
        het_len = random.randint(0, max_het)
        num_seqs = random.randint(1, max_seqs)
        seqs = [rand_seq(het_len + 1) for _ in range(num_seqs)]

        spacer_finder = FindSpacerCombo(10, NUM_THREADS, seqs, het_len)
        spacers = spacer_finder.run()

        spacer_seqs = get_rand_hetero_seqs(seqs, spacers)

        for i, spacer in enumerate(spacer_seqs):
            assert len(spacer) == spacers[i]

    seqs = [
        'ATGCTA',
        'CAGCTA',
        'ACTAGT',
        'ATGCTA',
        'TAGCTG'
    ]

    het_spacers = [0, 1, 2, 3, 4]

    spacer_seqs = get_rand_hetero_seqs(seqs, het_spacers)

    for i, spacer in enumerate(spacer_seqs):
        assert len(spacer) == het_spacers[i]


def test_simple_seq_generation() -> None:
    gen_runtime = 10

    f_binding = ['GAGTAGAACTGGGAGCTCATTCGA', 'TCGGTAAACGGCGCGCGTCTCCTA',
                 'AGTTCACTGGCCGAATACTCTCGT', 'CAGGCCTTAGGTGCCTGGCAGCCC']
    r_binding = ['GATCAGGCATCCCGGCGGCTTATA', 'ACTCCTCTATGTATACGGTACATG',
                 'GTACGCGTTAACAGAGGCCTTGGG', 'TAACGGTAACCGGTAATCTAATGT']
    f_adapters = ['AGAATCGCGGGAGCTAAGGTGGTA', 'TAGACCTACTAAGGTCGTGTATGG',
                  'GCGGTAGACACTGGCTTCGGGCGG', 'GGAACCATCGCAGATTTCGACAAT']
    r_adapters = ['ACCCCTAGCGCAATCCCCAACTGC', 'AATCAACCTACCGCGCGTACGTTC',
                  'ACGACGTACCGTAGCTCGCCTTCT', 'CTCGTCAGATCTCTAACTGGCGCA']
    f_spacers = [5, 0, 0, 2]
    r_spacers = [5, 1, 0, 2]
    oq = Queue()

    get_best_heterogeneity_spacer_seqs(f_adapters, f_binding,
                                       r_adapters, r_binding,
                                       f_spacers, r_spacers,
                                       gen_runtime, 5, oq)
    data = None
    while not oq.empty():
        data = oq.get()
        if isinstance(data, str):
            log.info(data)
        else:
            break
    best_pool = data
    assert isinstance(best_pool[1], PrimerPool)


def test_hetero_seq_generation_unthreaded() -> None:
    max_het = 12
    max_seqs = 5
    num_reps = 5
    max_adapter_len = 24
    max_seq_len = 24
    spacer_finder_runtime = 10
    seq_generation_runtime = 30

    for _ in range(num_reps):
        het_len = random.randint(0, max_het)
        num_seqs = random.randint(1, max_seqs)
        f_binding = [rand_seq(max_seq_len) for _ in range(num_seqs)]
        r_binding = [rand_seq(max_seq_len) for _ in range(num_seqs)]
        f_adapters = [rand_seq(max_adapter_len) for _ in range(num_seqs)]
        r_adapters = [rand_seq(max_adapter_len) for _ in range(num_seqs)]
        num_structs_to_view = num_seqs ** 2 // 2

        spacer_finder = FindSpacerCombo(spacer_finder_runtime, NUM_THREADS,
                                        f_binding, het_len)
        f_spacers = spacer_finder.run()
        spacer_finder = FindSpacerCombo(spacer_finder_runtime, NUM_THREADS,
                                        r_binding, het_len)
        r_spacers = spacer_finder.run()

        oq = Queue()

        get_best_heterogeneity_spacer_seqs(f_adapters, f_binding,
                                           r_adapters, r_binding,
                                           f_spacers, r_spacers,
                                           seq_generation_runtime,
                                           num_structs_to_view, oq)
        data = None
        while not oq.empty():
            data = oq.get()
            if isinstance(data, str):
                log.info(data)
            else:
                break
        best_pool = data
        assert isinstance(best_pool[1], PrimerPool)


def test_simple_seq_generation_threaded() -> None:
    gen_runtime = 60

    f_binding = ['GAGTAGAACTGGGAGCTCATTCGA', 'TCGGTAAACGGCGCGCGTCTCCTA',
                 'AGTTCACTGGCCGAATACTCTCGT', 'CAGGCCTTAGGTGCCTGGCAGCCC']
    r_binding = ['GATCAGGCATCCCGGCGGCTTATA', 'ACTCCTCTATGTATACGGTACATG',
                 'GTACGCGTTAACAGAGGCCTTGGG', 'TAACGGTAACCGGTAATCTAATGT']
    f_adapters = ['AGAATCGCGGGAGCTAAGGTGGTA', 'TAGACCTACTAAGGTCGTGTATGG',
                  'GCGGTAGACACTGGCTTCGGGCGG', 'GGAACCATCGCAGATTTCGACAAT']
    r_adapters = ['ACCCCTAGCGCAATCCCCAACTGC', 'AATCAACCTACCGCGCGTACGTTC',
                  'ACGACGTACCGTAGCTCGCCTTCT', 'CTCGTCAGATCTCTAACTGGCGCA']
    f_spacers = [5, 0, 0, 2]
    r_spacers = [5, 1, 0, 2]

    final_pool = get_best_heterogeneity_spacer_seqs_threadable(
        f_adapters, f_binding, r_adapters, r_binding, f_spacers, r_spacers,
        gen_runtime, 5, NUM_THREADS)
    print(final_pool)


def test_hetero_seq_generation_threaded() -> None:
    max_het = 12
    max_seqs = 20
    num_reps = 5
    max_adapter_len = 24
    max_seq_len = 24
    spacer_finder_runtime = 10
    seq_generation_runtime = 10

    for _ in range(num_reps):
        het_len = random.randint(0, max_het)
        num_seqs = random.randint(1, max_seqs)
        f_binding = [rand_seq(max_seq_len) for _ in range(num_seqs)]
        r_binding = [rand_seq(max_seq_len) for _ in range(num_seqs)]
        f_adapters = [rand_seq(max_adapter_len) for _ in range(num_seqs)]
        r_adapters = [rand_seq(max_adapter_len) for _ in range(num_seqs)]
        num_structs_to_view = num_seqs ** 2 // 2

        spacer_finder = FindSpacerCombo(spacer_finder_runtime, NUM_THREADS,
                                        f_binding, het_len)
        f_spacers = spacer_finder.run()
        spacer_finder = FindSpacerCombo(spacer_finder_runtime, NUM_THREADS,
                                        r_binding, het_len)
        r_spacers = spacer_finder.run()

        final_pool = get_best_heterogeneity_spacer_seqs_threadable(
            f_adapters, f_binding, r_adapters, r_binding, f_spacers, r_spacers,
            seq_generation_runtime, num_structs_to_view, NUM_THREADS)
        assert final_pool


if __name__ == '__main__':
    standard_pytest_run(__file__)
