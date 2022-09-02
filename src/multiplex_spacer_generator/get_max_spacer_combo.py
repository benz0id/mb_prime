from multiprocessing import Queue
import os
from typing import List

from src.multiplex_spacer_generator.binding_align import NumpyBindingAlign
from logging import Logger

NOTHING_FOUND_MSG = 'nothing_found'
LOG_EVERY = 10000
DEBUG_LOGGING = False


def get_max_spacer_combo(binding_align: NumpyBindingAlign, out_queue: Queue,
                         num_to_iter: int, max_score: int,
                         do_random_generation: bool,
                         max_spacer_len: int = -1, return_best: bool = False) -> None:
    """Will iterate over the next <num_to_iter> spacer combos in
    <start_binding_align> or until a combo is found with score equal to
    <max_score>, and will place that combo in <out_queue>."""
    minimum_size_threshold = max_spacer_len != -1

    if do_random_generation:
        num_to_iter = 10 ** 20
        random_gen = ' - RANDOM GENERATION ENABLED'
        binding_align.random_spacer_maintain_size()
    else:
        random_gen = ''

    process_header = 'Child Process ' + str(os.getpid()) + ': '

    best_spacers_so_far = []

    def return_if_valid(cur_binding_align: NumpyBindingAlign,
                        best_spacers_so_far_lst: List[NumpyBindingAlign]) -> bool:
        score = cur_binding_align.find_min_div()
        if DEBUG_LOGGING:
            out_queue.put(''.join(
                [
                    process_header, 'Evaluating ', str(cur_binding_align), '\n',
                    'Score: ', str(score), '\n',
                    'Max Len: ', str(max(binding_align.get_spacer_sizes()))
                ]
            ))
        if score == max_score:
            # We don't care about the max len, so return the combo.
            if not minimum_size_threshold:
                out_queue.put(cur_binding_align.get_spacer_sizes())
                return True
            # Is the combo under the max spacer len threshold? If not is it
            # closer than anything found so far?
            else:
                best_so_far_present = bool(best_spacers_so_far_lst)
                # Is it under the threshold?
                if max(binding_align.get_spacer_sizes()) <= max_spacer_len:
                    out_queue.put(cur_binding_align.get_spacer_sizes())
                    return True
                # Is it better than anything found so far?
                elif not best_so_far_present or \
                    max(cur_binding_align.get_spacer_sizes()) < \
                    max(best_spacers_so_far_lst[0].get_spacer_sizes()):
                    best_spacers_so_far.insert(0, cur_binding_align)
                    return True
        return False

    out_queue.put(''.join([
        process_header, 'Beginning iteration over ', str(num_to_iter),
        ' spacers, starting at ', str(binding_align.get_spacer_sizes()), '. ',
        random_gen
    ]))

    # Check initial align for maximal score.
    combo_found = return_if_valid(binding_align, best_spacers_so_far)

    # Check specified region for maximal score.
    i = 0
    for i in range(num_to_iter - 1):
        detailed_logging = i % LOG_EVERY == 0

        # Continue to next spacer combo that maintains the total size.
        try:
            if not do_random_generation:
                binding_align.incr_spacer_maintain_size()
            else:
                binding_align.random_spacer_maintain_size()
        except StopIteration:
            break

        # Do we only care about spacers with short lengths and is this spacer
        # over the max spacer length? If so, skip it.
        spacer_size_too_large = max(binding_align.get_spacer_sizes()) > \
                                max_spacer_len
        if not return_best and minimum_size_threshold and spacer_size_too_large:
            continue

        # If the score is equal to the maximal score, add it to the queue.
        combo_found = return_if_valid(binding_align, best_spacers_so_far)
        if detailed_logging:
            out_queue.put(''.join(
                [
                    process_header, 'Completed iteration over ', str(i),
                    ' combos. Current: ', str(binding_align.get_spacer_sizes()),
                    ' Score: ', str(binding_align.find_min_div()), '.'
                ]
            ))

    if not combo_found:
        out_queue.put(''.join(
            [
                process_header, 'Finished iteration without finding a valid '
                                'spacer combo.\n\tFinal combo checked: ',
                str(binding_align.get_spacer_sizes()), ' - ',
                str(binding_align.find_min_div()), '\n\t Num spacers checked: ',
                str(i + 1)
            ]))

        # Nothing with a maximal score was found, return nothing found message.
        out_queue.put(NOTHING_FOUND_MSG)

    # We didn't find anything matching our max length criteria, but return
    # the best one anyways.
    if return_best and combo_found and best_spacers_so_far:
        out_queue.put(best_spacers_so_far[0])

