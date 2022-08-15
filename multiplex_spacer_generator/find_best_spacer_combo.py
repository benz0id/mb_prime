import logging
from time import sleep, time
from typing import List, Tuple, Iterable
from multiprocessing import Process, Manager, Queue
from tqdm import tqdm
from multiplex_spacer_generator.binding_align import SpacerCombo, \
    compute_column_score, incr_repr, NumpyBindingAlign, smart_incr, \
    get_min_pos_repr, get_time_string
from multiplex_spacer_generator.get_max_spacer_combo import \
    get_max_spacer_combo, NOTHING_FOUND_MSG
import heapq

TICK_SPEED = 1000
MIN_PER_PROCESS = 1000
MAX_BEFORE_RANDOM = 12 ** 8
SHORT_TIME_THRESHOLD = 5 * 60


def test_combo(spacer_combo: SpacerCombo, seqs: list[str],
               num_hetero: int) -> None:
    """Tests whether the given combo has maximal diversity."""
    max_score = NumpyBindingAlign(seqs, [num_hetero] * len(seqs), num_hetero).\
        find_min_div()
    acc_score = NumpyBindingAlign(seqs, spacer_combo, num_hetero).\
        find_min_div()
    assert max_score == acc_score


def time_prog_bar(stop_time: int) -> None:
    """Shows a progress bar with the remaining time."""
    for _ in tqdm(range(int(time()), stop_time)):
        sleep(1)


log = logging.getLogger('root')


def delay(miliseconds: int) -> None:
    # Convert to seconds.
    sleep(miliseconds / 1000)


def tick_wait() -> None:
    delay(TICK_SPEED)


class FindSpacerCombo:
    """Responsible for finding the spacer combo with the smallest average length
    in the allotted time.

    === Private Attributes ===

    # Provided parameters.
    _max_threads: Maximum number of processes to run at any given time.
    _seqs: Binding regions, 5'-3'.
    _hetero_region_size: The length of the heterogeneity region.
        0 <= spacer_len <= _hetero_region_size for each spacer.
    _max_spacer_size: The maximum size of any spacer to be returned.
    _min_per_process: The minimum number of tasks that must be assigned to any
        given process

    # Derived parameters.
    _end_time: The time after which this class will return the best spacer combo
        it has found.
    _max_score: Maximum possible diversity score attainable by a spacer.
    _do_random_generation: When number of possible spacer combos becomes too
        large, threads are instructed to randomly select them rather than follow
        a defined set.


    # Runtime Storage.
    _best_combos: A heap of the best spacers found, sorted by total length.

    # Multiprocessing
    _threads: Processes spawned by this class.
    _thread_manager: Manages communication between spawned processes.
    _thread_out_queue: Receives information from the spawned processes.



    """

    # Provided parameters.
    _max_threads: int
    _seqs: List[str]
    _hetero_region_size: int
    _max_spacer_size: int
    _min_per_process: int

    # Derived parameters.
    _end_time: int  # Seconds since epoch.
    _max_score: int
    _do_random_generation: int

    # Runtime Storage.
    _best_combos: List[Tuple[int, SpacerCombo]]

    # Multiprocessing
    _threads: List[Process]
    _thread_manager: Manager
    _thread_out_queue: Queue

    def __init__(self, runtime: int, max_threads: int, seqs: Iterable[str],
                 hetero_region_size: int, max_spacer_size: int = -1,
                 min_per_process: int = MIN_PER_PROCESS) -> None:
        """Initialises this class using the given parameters.
        runtime: allowed runtime in seconds."""

        self._max_threads = max_threads
        self._seqs = [seq for seq in seqs]
        self._hetero_region_size = hetero_region_size
        self._min_per_process = min_per_process

        self._end_time = int(time()) + runtime
        self._max_score = compute_column_score([len(self._seqs), 0, 0, 0, 0],
                                               len(self._seqs))

        num_pos_spacers = hetero_region_size ** len(self._seqs)
        self._do_random_generation = num_pos_spacers > MAX_BEFORE_RANDOM
        if self._do_random_generation:
            log.info(''.join([
                'Number of possible spacers exceeds threshold for ordered '
                'iteration.',
                '\n\t', str(num_pos_spacers), ' > ', str(MAX_BEFORE_RANDOM)
            ]))
            log.warning('RANDOM GENERATION ENABLED - '
                        'SPACER QUALITY MAY BE LOWER')

        self._best_combos = []

        self._threads = []
        self._thread_manager = Manager()
        self._thread_out_queue = Manager().Queue()

    def run(self) -> SpacerCombo:
        """Runs for the specified amount of time, returning the best spacer
        combo found."""
        max_total_len = len(self._seqs) * self._hetero_region_size
        if self._end_time - time() < SHORT_TIME_THRESHOLD and \
                not self._do_random_generation:
            total_len = max_total_len
        else:
            total_len = (max_total_len // 3) * 2

        if self._do_random_generation:
            incr_size = len(self._seqs) // 5
        else:
            incr_size = len(self._seqs)

        # For every possible total spacer length, starting at the maximum,
        # decreasing by the number of seqs in the alignment each time.
        while 0 <= total_len <= max_total_len:
            log.info('    === Beginning Search for Spacer Combo at Length ' +
                     str(total_len) + ' ===    ')
            log.info('Generating Threads.')
            self._gen_threads(total_len)
            log.info('Starting Threads.')
            self._start_threads()

            # Couldn't find spacer combo.
            if not self._wait_for_result_or_time_expiry():
                total_len += incr_size // 3
            else:
                total_len -= incr_size
            if self._time_expired():
                log.info('Stopping search for spacer combos.')
                break
            else:
                log.info('Search for Spacer Combo at Length ' +
                str(total_len) + ' complete.\nTime Remaining: ' +
                         get_time_string(self._end_time - time()))

        if not self._best_combos:
            raise RuntimeError('Failed to find a primer set in the allotted '
                               'time. Try allotting more time or entering a '
                               'less complicated input.')

        log.info('Returning the best found combo found: ' +
                 str(self._best_combos[0]))
        return self._best_combos[0][1]

    def _wait_for_result_or_time_expiry(self) -> bool:
        """Waits for valid data from threads. Returns True iff a spacer has been
        found, returns false if the timer has expired or the child processes
        failed to find a valid set."""
        t0 = time()
        continuation_permitted = True
        num_threads = len(self._threads)
        num_failures = 0
        combo_found = False

        while continuation_permitted:
            tick_wait()

            # Handle output from child procs.
            while not self._thread_out_queue.empty():
                data = self._thread_out_queue.get()

                if data == NOTHING_FOUND_MSG:
                    num_failures += 1

                # Log data that's been piped out.
                elif isinstance(data, str):
                    log.info(data.strip())

                elif not combo_found and isinstance(data, list):
                    spacer_tup = sum(data), data[:]
                    log.info('Time before receiving first combo:'
                             + get_time_string(time() - t0))
                    log.info('Combo received: ' + str(spacer_tup))
                    heapq.heappush(self._best_combos, spacer_tup)
                    continuation_permitted = False
                    combo_found = True
                    test_combo(data, self._seqs, self._hetero_region_size)
                    break

                else:
                    print(data)

            if num_failures == num_threads:
                log.info('All threads failed to find a valid combo.')
                continuation_permitted = False

            if self._time_expired():
                log.info('Time limit reached to find a combo.')
                continuation_permitted = False

        # Kill all threads.
        log.info('Killing remaining threads.')
        for thread in self._threads:
            thread.terminate()

        while not self._thread_out_queue.empty():
            self._thread_out_queue.get()

        if combo_found:
            log.info(combo_str(spacer_tup[1], self._seqs))

        return combo_found

    def _time_expired(self) -> bool:
        """Returns whether this class has time remaining."""
        return self._end_time < time()

    def _show_time_prog_bar(self) -> None:
        """Starts a thread that shows a progress bar with the remaining time."""
        thread = Process(target=time_prog_bar,
                         args=tuple([self._end_time]))
        thread.start()

    def _start_threads(self) -> None:
        """Start the child processes started by this class."""
        for thread in self._threads:
            thread.start()

    def _gen_threads(self, total_len: int) -> None:
        """Generates the threads that will iterate over all spacer combos."""
        self._threads = []
        thread_str = ''.join(
            [
                '\t Thread #n : ', 'starting_spacer',
                ' - ', 'num_to_iter', '\n'
            ]
        )
        thread_num = 0
        # Convert each designated region into its own thread.
        for num_to_iter, starting_spacer in self._get_thread_params(total_len):

            thread_str += ''.join(
                [
                    '\t Thread #', str(thread_num), ': ', str(starting_spacer),
                    ' - ', str(num_to_iter), '\n'
                ]
            )
            binding_align = NumpyBindingAlign(self._seqs, starting_spacer,
                                              self._hetero_region_size)
            thread = Process(target=get_max_spacer_combo,
                             args=(binding_align, self._thread_out_queue,
                                   num_to_iter, self._max_score,
                                   self._do_random_generation))
            self._threads.append(thread)
            thread_num += 1

        log.info(thread_str)

    def _get_num_pos(self, total_len: int) -> int:
        """Returns the number possible of spacer alignments that contain a total
        of <total_bases> in their heterogeneity region."""
        repr = get_min_pos_repr(total_len, len(self._seqs),
                                self._hetero_region_size)

        count = 1
        while smart_incr(repr, self._hetero_region_size):
            count += 1

        return count

    def _get_nth_pos(self, nths: List[int], total_len: int) -> List[
        SpacerCombo]:
        """Returns the number possible of spacer alignments that contain a total
        of <total_bases> in their heterogeneity region."""
        nths = nths[:]

        repr = get_min_pos_repr(total_len, len(self._seqs),
                                self._hetero_region_size)

        nths.sort(reverse=True)
        nth_spacer_combos = []
        count = 0
        if nths and nths[-1] == count:
            nth_spacer_combos.append(repr[:])
            nths.pop()

        # Fetch each of the nth specified values.
        while smart_incr(repr, self._hetero_region_size):
            count += 1
            if nths and nths[-1] == count:
                nth_spacer_combos.append(repr[:])
                nths.pop()

        if len(nths) > 0:
            raise ValueError('Inappropriate inputs received.')

        return nth_spacer_combos

    def _get_thread_params(self, total_len: int) \
            -> List[Tuple[int, SpacerCombo]]:
        """Returns a list of initialisation parameters for (at most)
        <self._max_threads>. Each will specify the starting spacer combo and how
        many more to iterate beyond that first spacer combo. If all threads run
        to completion, every heterogeneity spacer with <total_len> will have
        been iterated over."""
        t0 = time()
        # If spacers are to be generated randomly, we needn't delegate regions.
        if self._do_random_generation:
            return [(-1, get_min_pos_repr(total_len, len(self._seqs),
                                          self._hetero_region_size))
                    for _ in range(self._max_threads)]

        num_pos = self._get_num_pos(total_len)
        if num_pos == 0:
            raise ValueError('There exist no possible spacer combos that '
                             'produce the specified length.')

        # Try to use max number of threads. Decrease to avoid excessive process
        # spawning.
        num_threads = self._max_threads
        while num_pos // num_threads < self._min_per_process \
                and not num_threads == 1:
            num_threads -= 1
        base_responsibility = num_pos // num_threads

        # Assign responsibility, account for task remainder.
        thread_responsibilities = [base_responsibility] * num_threads
        thread_responsibilities[-1] += num_pos % num_threads

        start_points = [sum(thread_responsibilities[0:i])
                        for i in range(0, len(thread_responsibilities))]

        # Each thread i will perform thread_responsibilities[i] tasks starting
        # at and including start_points[i].
        start_combos = self._get_nth_pos(start_points, total_len)

        # Pair results into tuples.
        rtrn = [tuple((thread_responsibilities[i], start_combos[i]))
                for i in range(len(start_combos))]

        td = time() - t0
        log.info('Time spent getting thread params:' + get_time_string(td))

        return rtrn


def combo_str(spacers: SpacerCombo, seqs: List[str]) -> str:
    """Coverts the spacer combo into a string."""
    s = str(spacers) + '\n'
    for i in range(len(spacers)):
        s += ' '.join('-' * spacers[i] + seqs[i] + '\n')

    return s