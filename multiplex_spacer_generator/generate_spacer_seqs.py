import os
from typing import List, Any, Callable
import math
import logging
from multiplex_spacer_generator.binding_align import SpacerCombo
import random
from multiprocessing import Queue, Manager, Process
from time import time, sleep
from statistics import mean
from multiplex_spacer_generator.primer_pool import PrimerPool, get_all_dgs

log = logging.getLogger('root')

BASE_TO_INT = {
    '-': 0,
    'A': 1,
    'T': 2,
    'G': 3,
    'C': 4
}
INT_TO_BASE = {v: k for k, v in BASE_TO_INT.items()}

NUM_ITER_BEFORE_CHECK = 100

ALERT_EVERY = 1000

MIN_DG = -10000000000000000000


def gen_heterogeneity_alignment(seqs: List[str], spacers: SpacerCombo) \
        -> List[List[str]]:
    """Generates an empty heterogeneity alignment, with '-' representing
    unassigned bases that lie in a heterogeneity spacer."""
    het_len = max(spacers)
    seq_len = len(seqs)
    arr = [['-'] * het_len for _ in range(seq_len)]

    for i, seq in enumerate(seqs):
        spacer_len = spacers[i]
        arr[i][spacer_len: het_len] = seq[:het_len - spacer_len]

    return arr


def get_allowable_bases(aln: List[List[str]], column: int) -> List[int]:
    """Returns the number of each base that could be inserted while maintaining
    heterogeneity."""
    next_4 = math.ceil(len(aln) / 4)
    allowable = [next_4] * 4

    for row in range(len(aln)):
        base = aln[row][column]
        if base == '-':
            continue
        ind = BASE_TO_INT[base] - 1
        allowable[ind] -= 1

    for num_base in allowable:
        assert num_base >= 0

    return allowable


def sel_random_base(allowable_bases: List[int]) -> str:
    """Selects and returns a base from the given allowable bases."""
    base_inds = []
    for ind, base in enumerate(allowable_bases):
        if base:
            base_inds.append(ind)
    rand_base_ind = random.choice(base_inds)
    allowable_bases[rand_base_ind] -= 1
    return INT_TO_BASE[rand_base_ind + 1]


def get_rand_hetero_seqs(seqs: List[str], spacers: SpacerCombo) -> List[str]:
    """Generates heterogeneity of length <spacers> for the given <seqs>."""
    aln = gen_heterogeneity_alignment(seqs, spacers)

    for col in range(len(aln[0])):
        allowable_bases = get_allowable_bases(aln, col)

        for row in range(len(aln)):
            if aln[row][col] in 'ATCG':
                continue
            aln[row][col] = sel_random_base(allowable_bases)

    seqs = [''.join(aln[ind][:spacers[ind]]) for ind in range(len(spacers))]

    return seqs


def get_new_primer_pool(f_5p: List[str], f_binding: List[str],
                        r_5p: List[str], r_binding: List[str],
                        f_spacers: SpacerCombo, r_spacers: SpacerCombo) \
        -> PrimerPool:
    """Returns a primer pool with randomly generated heterogeneity spacers."""
    f_spacer_seqs = get_rand_hetero_seqs(f_binding, f_spacers)
    r_spacer_seqs = get_rand_hetero_seqs(r_binding, r_spacers)
    return PrimerPool(f_5p, f_spacer_seqs, f_binding, r_5p, r_spacer_seqs,
                      r_binding)


def get_best_heterogeneity_spacer_seqs(
        f_5p: List[str], f_binding: List[str], r_5p: List[str],
        r_binding: List[str], f_spacers: SpacerCombo, r_spacers: SpacerCombo,
        allowed_seconds: int, num_structs_to_avg: int,
        out_queue: Queue) -> None:
    """Runs for <allowed_seconds> after which the best PrimerPool will be
    returned."""
    best_set = None
    highest_avg_dg = MIN_DG
    num_iter = 0

    process_header = 'Child Process ' + str(os.getpid()) + ': '

    end_time = time() + allowed_seconds

    while end_time > time():
        for _ in range(NUM_ITER_BEFORE_CHECK):

            if num_iter % ALERT_EVERY == 0:
                msg = ''.join(
                    [
                        process_header, str(num_iter), ' completed. Best dg:',
                        str(int(highest_avg_dg))
                    ]
                )
                out_queue.put(msg)

            seq_pool = get_new_primer_pool(f_5p, f_binding, r_5p, r_binding,
                                           f_spacers, r_spacers)
            seqs = [str(seq) for seq in seq_pool.get_all_seqs()]
            struct_dgs = get_all_dgs(seqs)

            mean_dg = mean(struct_dgs[-num_structs_to_avg:])

            if mean_dg > highest_avg_dg:
                msg = ''.join(
                    [
                        process_header, ' New best set found. [',
                        str(int(highest_avg_dg)), ' --> ', str(int(mean_dg)),
                        ']'
                    ]
                )
                highest_avg_dg = mean_dg
                best_set = seq_pool
                out_queue.put(msg)
            num_iter += 1
    msg = ''.join(
        [
            process_header, ' Returning best set found - ',
            str(int(highest_avg_dg)),
        ]
    )
    out_queue.put(msg)
    out_queue.put((highest_avg_dg, best_set))


def get_best_heterogeneity_spacer_seqs_threadable(
        f_5p: List[str], f_binding: List[str], r_5p: List[str],
        r_binding: List[str], f_spacers: SpacerCombo, r_spacers: SpacerCombo,
        allowed_seconds: int, num_structs_to_avg: int,
        num_threads: int) -> PrimerPool:
    """Returns the primer pool least predisposed to forming very stable dimer
    structures."""
    manager = Manager()
    thread_out_queue = manager.Queue()
    fin_time = time() + allowed_seconds

    log.info('   ==== Beginning Search for Least Dimer-Prone Heterogeneity '
             'Spacer Sequences ====')

    threads = []
    for _ in range(num_threads):
        threads.append(Process(target=get_best_heterogeneity_spacer_seqs,
                               args=(f_5p, f_binding, r_5p, r_binding,
                                     f_spacers, r_spacers, allowed_seconds,
                                     num_structs_to_avg, thread_out_queue)))

    for thread in threads:
        thread.start()

    def all_joined() -> bool:
        for thread in threads:
            if thread.is_alive():
                return False
        return True

    time_left = True
    all_procs_joined = False
    all_data_received = False
    num_data_rec = 0

    # Record output data.
    rslts = []
    while time_left or not all_procs_joined or not all_data_received:
        while not thread_out_queue.empty():
            data = thread_out_queue.get()
            if isinstance(data, str):
                log.info(data)
            else:
                num_data_rec += 1
                rslts.append(data)
        time_left = fin_time > time()
        all_procs_joined = all_joined()
        all_data_received = num_data_rec == num_threads


        sleep(0.1)


    max_set = max(rslts, key=lambda a: a[0])[1]

    return max_set





