from multiprocessing import Queue
import os
from multiplex_spacer_generator.binding_align import NumpyBindingAlign
from logging import Logger

NOTHING_FOUND_MSG = 'nothing_found'


def get_max_spacer_combo(binding_align: NumpyBindingAlign, out_queue: Queue,
                         num_to_iter: int, max_score: int) -> None:
    """Will iterate over the next <num_to_iter> spacer combos in
    <start_binding_align> or until a combo is found with score equal to
    <max_score>, and will place that combo in <out_queue>."""
    process_header = 'Child Process ' + str(os.getpid()) + ': '
    combo_found = False

    out_queue.put(''.join([
        process_header, 'Beginning iteration over ', str(num_to_iter),
        ' spacers, starting at ', str(binding_align.get_spacer_sizes()), '.'
    ]))

    # Check initial align for maximal score.
    score = binding_align.find_min_div()
    if score == max_score:
        out_queue.put(binding_align.get_spacer_sizes())
        combo_found = True

    # Check specified region for maximal score.
    for _ in range(num_to_iter - 1):
        # Continue to next spacer combo that maintains the total size.
        binding_align.incr_spacer_maintain_size()
        # If the score is equal to the maximal score, add it to the queue.
        score = binding_align.find_min_div()
        if score == max_score:
            out_queue.put(binding_align.get_spacer_sizes())
            combo_found = True
    if not combo_found:
        out_queue.put(''.join(
            [
                process_header, 'Finished iteration without finding a valid '
                                'spacer combo.\n\tFinal combo checked: ',
                str(binding_align.get_spacer_sizes()), ' - ',
                str(binding_align.find_min_div())
            ]))

        # Nothing with a maximal score was found, return nothing found message.
        out_queue.put(NOTHING_FOUND_MSG)

