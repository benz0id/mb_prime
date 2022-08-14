import os
import timeit
from pathlib import Path
from sys import stdout
from typing import List, Tuple, Union
import math
import numpy as np
from copy import deepcopy
from random import randint, choices, choice
from statistics import stdev
import time
from multiprocessing import Process, Manager, Queue
import logging as lg

logger = lg.RootLogger(level=0)
logger.addHandler(lg.StreamHandler(stdout))

# Apple trees and goat cheese.
BASE_TO_INT = {
    '-': 0,
    'A': 1,
    'T': 2,
    'G': 3,
    'C': 4
}


def get_min_pos_repr(total: int, len: int, B: int) -> List[int]:
    """Returns a list containing the minimal radix <B> + 1 representation with
    sum == total"""
    lst = [0] * len
    i = 0
    while total != 0:
        amt = min(total, B)
        lst[i] = amt
        total -= amt
        i += 1
    return lst


def incr_repr(repr: List[int], B: int) -> bool:
    """Increments a representation <repr> interpreted to be radix <B> + 1 with
    the least significant digit at <repr>[0]."""
    for i, v in enumerate(repr):
        if v < B:
            repr[i] += 1
            return True
        else:
            repr[i] = 0
            if i == len(repr) - 1:
                return False


def smart_incr(lst: List[int], max_val: int) -> bool:
    """Increments radix <max_val> + 1 representaion <lst> S.T. sum(lst) remains
    unnchanged."""
    def getmin() -> int:
        for i in range(len(lst) - 1):
            if lst[i] > 0 and lst[i + 1] < max_val:
                return i
        return -1

    min_ind = getmin()

    if min_ind == -1:
        return False

    lst[min_ind] -= 1
    lst[min_ind + 1] += 1

    sums = sum(lst[:min_ind + 1])
    for i in range(min_ind + 1):
        lst[i] = 0

    i = 0
    while sums != 0:
        amt = min(sums, max_val)
        lst[i] = amt
        sums -= amt
        i += 1

    return True


def compute_column_score(bases: Union[List[int], np.ndarray],
                         total: int) -> int:
    """Given the base frequency of each type of base in an alignment column,
    ordered: unassigned A, T, C, G, (matching indexing with BASE_TO_INT),
    returns some score quantifying the level (or potential level) of diversity
    in that column. Total is the total number of bases in that column.

    Precondition:
        <bases> must match the formatting given, being of length 5."""

    b_inds = [0, 1, 2, 3]
    b_vals = np.array(bases[1:], dtype=np.int8)

    # Sort from least to most occurrences.
    list.sort(b_inds, key=lambda a: b_vals[a])

    unassigned = bases[0]
    while unassigned > 0:
        # Assign base to the lowest abundance base.
        b_vals[b_inds[0]] += 1
        unassigned -= 1

        # Adjust position of incremented base to make sure indices are properly
        # sorted.
        for i in range(3):
            if b_vals[b_inds[i]] > b_vals[b_inds[i + 1]]:
                b_inds[i], b_inds[i + 1] = b_inds[i + 1], b_inds[i]
            else:
                break
    # Calculate percent representation.
    percent_reps = [max(val / total * 100, 1) for val in b_vals]

    return int(np.prod(percent_reps))


SpacerCombo = List[int]


class NumpyBindingAlign:
    """A class designed to enable rapid comparison of diversity at within a
    short alignment.

    _og_align:
        The original alignment of DNA sequences.

    _numpy_align:
        An 2D array representation of the alignment, with bases being converted
        to numbers.

    _spacer_sizes:
        The sizes of the heterogeneity spacers used for each of the sequences in
        the alignment. This should be an alias of some other list.

    _num_seqs:
        The number of sequences in the alignment.

    _num_heterogeneity:
        Size of the heterogeneity region."""

    _og_align: List[str]

    _numpy_align: np.ndarray

    _spacer_sizes: List[int]

    _num_seqs: int

    _num_hetero: int

    def __init__(self, align: List[str], spacer_sizes: List[int],
                 num_hetero: int) -> None:
        """Converts the alignment to a numpy array and stores it internally.

        Preconditions:
            Number of bases of heterogeneity must be less than length of the
            shortest sequence in alignment."""
        self._og_align = deepcopy(align)
        self._num_hetero = num_hetero
        self._num_seqs = len(align)

        lens = [len(seq) for seq in align]
        min_len = min(lens)
        max_len = max(lens)

        if min_len < num_hetero:
            raise ValueError("Number of bases of heterogeneity must be less "
                             "than the length of the shortest sequence in "
                             "alignment.")

        # Create integer representation of the alignment.
        int_rep_arr = []
        for i, seq in enumerate(align):
            int_rep_arr.append([])
            j = max_len
            for j, char in enumerate(seq):
                num_rep = BASE_TO_INT[char]
                int_rep_arr[i].append(num_rep)
            for k in range(j, max_len - 1):
                int_rep_arr[i].append(BASE_TO_INT['-'])

        # Convert to numpy array.
        self._numpy_align = np.array(int_rep_arr, dtype=np.int8)

        # Create spacer sizes. ALIASING!
        if len(spacer_sizes) != self._num_seqs:
            raise ValueError('Spacer list should be same length as alignment.')

        self._spacer_sizes = spacer_sizes

    def get_spacer_sizes(self) -> SpacerCombo:
        """Returns a copy of the spacer combo currently stored by this binding
        align."""
        return self._spacer_sizes[:]

    def find_min_div(self) -> int:
        """Returns the lowest diversity value found in the  first
        <self._num_hetero> bases."""
        assert min(self._spacer_sizes) >= 0
        assert max(self._spacer_sizes) <= self._num_hetero

        # Set to theoretical max.
        min_score = 25 ** 4
        # Find diversity score of each column.
        for col in range(self._num_hetero):
            # Number of each type of base in the current column,
            #                     A, T, C, G, unassigned
            # In accordance with BASE_TO_INT
            num_bases = np.array([0, 0, 0, 0, 0], dtype=np.int32)

            # Count number of each nucleotide in each sequence within the
            # column.
            for b in range(self._num_seqs):
                # Base lies in heterogeneity spacer.
                if col < self._spacer_sizes[b]:
                    num_bases[0] += 1
                    continue

                # Accounting for the shift introduced by heterogeneity spacers.
                shifted_col = col - self._spacer_sizes[b]
                # Increment base count.
                num_bases[self._numpy_align[b, shifted_col]] += 1

            # Compute score.
            score = compute_column_score(num_bases, self._num_seqs)
            if score < min_score:
                min_score = score

        return min_score

    def incr_spacer_maintain_size(self) -> None:
        """Increments the spacer while maintaining the total length of spacers.
        """
        if not smart_incr(self._spacer_sizes, self._num_hetero):
            raise StopIteration

    def random_spacer_maintain_size(self) -> None:
        """Sets the spacer to another random spacer with the same size as the
        current one."""
        total_size_to_add = sum(self._spacer_sizes)
        self._spacer_sizes = [0] * len(self._spacer_sizes)
        inds = list(range(len(self._spacer_sizes)))

        while total_size_to_add > 0:
            ind = choice(inds)
            self._spacer_sizes[ind] += 1
            if self._spacer_sizes[ind] == self._num_hetero:
                inds.remove(ind)
            total_size_to_add -= 1

    def get_mean_spacer_size(self) -> float:
        """Returns the mean spacer size"""
        return sum(self._spacer_sizes) / len(self._spacer_sizes)

    def get_spacer_std(self) -> float:
        """Returns the deviation amoung the spacers."""
        return stdev(self._spacer_sizes)

    def incr_spacer_sizes(self, num: int = 1) -> None:
        """Increments the spacer sizes."""
        for _ in range(num):
            incr_repr(self._spacer_sizes, self._num_hetero)


def gen_random_seq_str(length: int, GC: float = 0.5) -> str:
    """Generates a random Seq of <length>"""
    AT = 1 - GC
    bases = ['A', 'T', 'C', 'G']

    return ''.join(choices(bases, weights=[AT / 2, AT / 2, GC / 2, GC / 2],
                           k=length))

def get_runtime(num_seqs: int, num_hetero: int, silent: bool = False) -> float:
    """Gets the runtime for some number of find_max_div operations. Returns the
    average runtime of each operation."""

    # Inputs for a standard use case.
    len_seqs = 20
    reps = 1000
    seqs = [gen_random_seq_str(len_seqs) for _ in range(num_seqs)]
    spacer_sizes = [0 for _ in range(num_seqs)]

    aln = NumpyBindingAlign(seqs, spacer_sizes, num_hetero)

    def run() -> None:
        for i in range(len(spacer_sizes)):
            spacer_sizes[i] = randint(0, num_hetero)
        score = aln.find_min_div()

    t = timeit.timeit(lambda: run(), number=reps)
    sec_per_op = t / reps
    op_per_sec = 1 / sec_per_op
    s = ['Speed test runtime:', str(t), '\n',
         'Sec per operation:', str(sec_per_op), '\n',
         'Operations per sec:', str(op_per_sec)
         ]
    if not silent:
        logger.info(' '.join(s))
    return sec_per_op

def get_time_string(secs: int or float) -> str:
    """Convert the seconds value into a more readable time string."""
    h = secs // (60 ** 2)
    secs -= h * 60 ** 2
    m = secs // 60
    secs -= m * 60
    secs = int(secs)
    s = []
    if h > 0:
        s += [str(h), 'hours,']
    if m > 0:
        s += [str(m), 'minutes, and', str(secs), 'seconds']
    elif secs > 0:
        s += [str(secs), 'seconds']
    else:
        s += ['less than a second']

    return ' '.join(s)


def get_num_comp(num_seqs: int, len_hetero: int, num_reps: int = 1,
                 silent: bool = False,  num_cores: int = 1) -> float:
    """Prints the time required to fully explore all possible heterogeneity
    spacers as described"""
    time_per_comp = get_runtime(num_seqs, len_hetero)

    t = (len_hetero + 1) ** num_seqs * time_per_comp

    ms = ''
    if num_cores > 1:
        t = t / num_cores
        ms = 'using ' + str(num_cores) + ' cores'

    s = [
        'To fully explore an alignment of ', str(num_seqs), ' sequences with a ',
        'heterogeneity region of size ', str(len_hetero), ' it would take ',
        get_time_string(t), ms,  '. '
    ]

    if num_reps > 1:
        s.extend([
            'Total time for all ', str(num_reps), ' replicates is ',
            get_time_string(num_reps * t), '.'
        ])
    if not silent:
        logger.info(''.join(s))
    return t * num_reps


def get_results(aln: NumpyBindingAlign, num: int,
                out_queue: Queue, time_queue: Queue) -> None:
    """Finds the worst spacer scores in the next <num> possible alignments of
     <aln> and stores them in <worst_scores>."""

    # Initialise for storing output data.
    worst_scores = np.empty(shape=[num], dtype=[
        ('mean_length', float), ('standard_deviation', float),
        ('worst_score', int)
    ])

    # Find worst scores.
    for i in range(num):
        me = aln.get_mean_spacer_size()
        sd = aln.get_spacer_std()
        score = aln.find_min_div()
        worst_scores[i] = (me, sd, score)

        aln.incr_spacer_sizes()

    pt0 = time.time()
    out_queue.put(worst_scores)
    pt = time.time() - pt0
    time_queue.put(pt)
    return


def divide_responsibility(num_tasks: int, num_procs: int) -> List[int]:
    """Divides the number of total tasks between each process as equally as
    possible.

    Precondition:
        num_tasks must be greater than num_procs"""
    if num_procs == 1:
        return [num_tasks]

    lst = []
    resp_per_proc = num_tasks // num_procs
    remainder = num_tasks % num_procs

    for _ in range(num_procs):
        lst.append(resp_per_proc)
    for i in range(remainder):
        lst[i] += 1

    return lst


def gen_csv(num_seqs: int, num_hetero: int, num_reps: int,
            output_filepath: Path, num_procs: int = 3) -> None:
    """Outputs several files containing the lowest diversity scores from
    randomly generated alignments matching the above specifications."""
    len_seqs = num_hetero
    if not output_filepath.exists():
        os.mkdir(output_filepath)

    ma = Manager()
    results_queue = ma.Queue()
    time_queue = ma.Queue()
    t0 = time.time()
    queue_put_time = 0
    mt = 0
    pt = 0

    print('Starting generation...')

    # For each replicate,
    for r in range(num_reps):

        # Create the csv to store the output data.
        out = output_filepath / ('R' + str(r) + '.csv')
        file = open(out, 'w')

        # Randomly generate some sequences to occupy the alignment and store
        # them.
        seqs = [gen_random_seq_str(len_seqs) for _ in range(num_seqs)]
        with open(output_filepath / ('R' + str(r) + '.txt'), 'w') as txt:
            for seq in seqs:
                txt.write(seq + '\n')

        # Create child processes and divide responsibility.
        spacer_sizes = [0 for _ in range(num_seqs)]
        aln = NumpyBindingAlign(seqs, spacer_sizes, num_hetero)
        total_tasks = (num_hetero + 1) ** num_seqs
        resps = divide_responsibility(total_tasks, num_procs)
        procs = []

        for tasks in resps:
            aln_copy = deepcopy(aln)
            procs.append(Process(target=get_results,
                                 args=(aln_copy, tasks, results_queue,
                                       time_queue)))
            # Increment the spacer sizes.
            aln.incr_spacer_sizes(num=tasks)

        pt0 = time.time()
        for proc in procs:
            proc.start()
        for proc in procs:
            proc.join()
        pt += time.time() - pt0

        for i in range(num_procs):
            queue_put_time += time_queue.get()

        mt0 = time.time()
        while not results_queue.empty():
            results = results_queue.get()
            for result in results:
                me = result['mean_length']
                sd = result['standard_deviation']
                score = result['worst_score']
                file.write(', '.join([str(me), str(sd), str(score)]) + '\n')
        file.close()
        mt += time.time() - mt0

    tf = time.time()
    time_string = ''.join([
        'Time spend adding to Queues: ', get_time_string(queue_put_time), '\n',
        'Time extracting from Queues and writing to files: ',
        get_time_string(mt), '\n',
        'Time spend in processes: ', get_time_string(pt), '\n',

    ])
    print('Complete. Runtime:', get_time_string(tf - t0), '\n' + time_string,
          '\n')


def org_set():
    NUM_SEQS = 5
    NUM_HETERO = 12
    NUM_REPS = 4
    PATH = Path('C:\\Users\\bfern\\OneDrive - University of Toronto\\Desktop\\Hetero Output\\Het12 Seq5')

    tpc = get_runtime(NUM_SEQS, NUM_HETERO)
    get_num_comp(NUM_SEQS, NUM_HETERO)
    gen_csv(NUM_SEQS, NUM_HETERO, NUM_REPS, PATH)

    NUM_SEQS = 6
    NUM_HETERO = 12
    NUM_REPS = 3
    PATH = Path('C:\\Users\\bfern\\OneDrive - University of Toronto\\Desktop\\Hetero Output\\Het12 Seq6')

    tpc = get_runtime(NUM_SEQS, NUM_HETERO)
    get_num_comp(NUM_SEQS, NUM_HETERO)
    gen_csv(NUM_SEQS, NUM_HETERO, NUM_REPS, PATH)

    NUM_SEQS = 7
    NUM_HETERO = 11
    NUM_REPS = 1
    PATH = Path('C:\\Users\\bfern\\OneDrive - University of Toronto\\Desktop\\Hetero Output\\Het11 Seq7')

    tpc = get_runtime(NUM_SEQS, NUM_HETERO)
    get_num_comp(NUM_SEQS, NUM_HETERO)
    gen_csv(NUM_SEQS, NUM_HETERO, NUM_REPS, PATH)



if __name__ == '__main__':
    NUM_PROCS = 8
    NUM_HETERO = 12
    NUM_REPS = 3
    NUM_SEQS = 9
    O_PATH = Path('C:\\Users\\bfern\\OneDrive - University of Toronto\\Desktop\\Hetero Output')

    # Runtime Estimate
    rt = 0
    rt += get_num_comp(NUM_SEQS, NUM_HETERO, NUM_REPS, silent=False,
                 num_cores=NUM_PROCS)
    print("Est. total runtime:", get_time_string(rt))

    #for NUM_SEQS in range(4, 9):
    #    PATH = O_PATH / ('Het5 Seq' + str(NUM_SEQS))
    #    get_num_comp(NUM_SEQS, NUM_HETERO, NUM_REPS)
    #    gen_csv(NUM_SEQS, NUM_HETERO, NUM_REPS, PATH, num_procs=NUM_PROCS)




