from hetero_spacer_generator.sequence_tools import is_degen
from presenters import ConsolePresenter
from demo_tools import *
from time import time

PRESENTER = ConsolePresenter()
hg = HeteroGen(presenter=PRESENTER)
input("Press enter to start.")


incomplete_forward_primer, forward_spacer = get_primer_and_spacers(hg,
                                                                   'forward',
                                                                   one_step=True)
incomplete_reverse_primer, reverse_spacer = get_primer_and_spacers(hg,
                                                                   'reverse',
                                                                   one_step=True)

while True:

    degen = False
    for primer in (incomplete_forward_primer, incomplete_reverse_primer):
        for seq in primer:
            if is_degen(seq):
                degen = True
                break


    hg.set_pairwise()

    rigour = int(while_not_valid(
        "Enter the rigour with which the program "
        "should search for possible combinations of spacers. "
        "Note: increasing the rigour leads to "
        "an exponential increase in runtime. Enter a value between 1 and 10:",
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






