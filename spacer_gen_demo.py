from hetero_spacer_generator import HeteroGen
from presenters import ConsolePresenter
from demo_tools import *
from threading import Thread
from time import time

PRESENTER = ConsolePresenter()
hg = HeteroGen(presenter=PRESENTER)
input("Press enter to start.")
succ = False

incomplete_forward_primer = MBPrimerBuilder()
forward_spacers = []

while not succ:
    print("Please enter below the components of your forward primer.\n"
          "The adapter and binding regions are mandatory.\n")

    incomplete_forward_primer = get_incomplete_primer()

    print("Please input the desired parameters for the construction of your "
          "forward primers.")

    spacer_length, num_hetero = get_params()
    hg.set_params(spacer_length, num_hetero)

    forward_spacers = hg.get_all_spacer_combos(
        incomplete_forward_primer.get_binding_seq())

    if not forward_spacers:
        print("Failed to find a way in which to align the given sequences. Try"
              " setting the max spacer length to a greater value.")
    else:
        succ = True

print("Spacers found. Press enter to select desired spacer combo.")

num_to_spacer = hg.visualise_spacer_alignments(forward_spacers,
                                               incomplete_forward_primer.
                                               get_binding_seq())

max_val = max(num_to_spacer.keys())

spacer_ind = while_not_valid("Enter the number of the primer alignment you'd "
                             "like to use:",
                             "Invalid input: Ensure the number you've entered "
                             "represents one of the primers displayed.\n "
                             "Enter the number of the primer alignment you'd "
                             "like to use:", RANGE, start=1, end=max_val)

forward_spacer = num_to_spacer[int(spacer_ind)]

# Please forgive heinous copy pasta coding - in a rush

succ = False
incomplete_reverse_primer = MBPrimerBuilder()
reverse_spacers = []

while not succ:
    print("Please enter below the components of your reverse primer.\n"
          "The adapter and binding regions are mandatory.\n")

    incomplete_reverse_primer = get_incomplete_primer()

    print("Please input the desired parameters for the construction of your "
          "reverse primers.")

    spacer_length, num_hetero = get_params()
    hg.set_params(spacer_length, num_hetero)

    reverse_spacers = hg.get_all_spacer_combos(
        incomplete_reverse_primer.get_binding_seq())

    if not reverse_spacers:
        print("Failed to find a way in which to align the given sequences. Try"
              " setting the max spacer length to a greater value.")
    else:
        succ = True

print("Spacers found. Press enter to select desired spacer combo.")

num_to_spacer = hg.visualise_spacer_alignments(reverse_spacers,
                                               incomplete_reverse_primer.
                                               get_binding_seq())

max_val = max(num_to_spacer.keys())

spacer_ind = while_not_valid("Enter the number of the primer alignment you'd "
                             "like to use:",
                             "Invalid input: Ensure the number you've entered "
                             "represents one of the primers displayed.\n "
                             "Enter the number of the primer alignment you'd "
                             "like to use:", RANGE, start=1, end=max_val)

reverse_spacer = num_to_spacer[int(spacer_ind)]
while True:

    rigour = int(while_not_valid("Enter the rigour with which the program should "
                                 "search for possible combinations of spacers. "
                                 "Note: increasing the rigour leads to "
                                 "an exponential increase in runtime. Enter a "
                                 "value between 1 and 10:",
                                 "Please enter a number between 1 and 10:",
                                 RANGE, start=1, end=10))

    hg.set_rigour(rigour)

    number_to_return = int(
        while_not_valid("Enter the number of sets of primer sets you'd like to "
                        "see.", "Please enter a number between 1 and 10:",
                                 RANGE, start=1, end=10))

    print("Generating primers. This may take some time...")

    t0 = time()

    psets = hg.get_hetero_seqs(incomplete_forward_primer, incomplete_reverse_primer,
                       forward_spacer, reverse_spacer, number_to_return)

    runtime = time() - t0

    input("Completed search in {time:.2f} seconds. Press enter t"
          "o see results.".format(time=runtime))

    for i in range(len(psets)):
        print("Set #", i + 1, ':')
        print(str(psets[i]))

    input("Press enter to generate a new set of primers.")






