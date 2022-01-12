from presenters import ConsolePresenter
from time import time
from hetero_spacer_generator.primer_tools import MBPrimerBuilder
from typing import Dict
import random
from hetero_spacer_generator.spacer_generator.hetero_spacer_generator import \
    HeteroGen
from pathlib import Path
from template_sequences import *

DESKTOP = Path('C:\\Users\\bfern\\OneDrive - University of Toronto\\Desktop')

# Defaults for generation of primers.
RIGOUR = 1
# The number of sets generated for binding sequence - adapter pair.
NUM_SETS = 1
# Output file for primers.
OUT_FILE = DESKTOP / Path(
    "Analysis Datasets\\Optimal Primers - Rigour\\Primer Sets")
# Print additional information
VERBOSE = True
# Format the output in fasta format. Required for TF analysis.
DO_FASTA = True
# Starting adapters and binding sequences. See available in
# meta_tools.template_sequences.py
GEN_SET = OPTIMAL_SET
# Whether output files should actually be made.
GEN_OUTPUT = True
# Whether you'd like seperate files for each binding seq - adapter pair.
MULTI_FILE = False


# I've left the above variables as default values in get_primer_sets, in case
# you'd like to vary them, as shown below.

# Put whatever you'd like to do in here.
def main():
    for rigour in range(-20, 20):
        filename = "Rigour " + str(rigour)
        print(filename)
        get_primer_sets(GEN_SET, NUM_SETS, filename=filename,
                        rigour=rigour)


def get_n_primer_sets(for_adapter: Seq, rev_adapter: Seq, for_binding: Seq,
                      rev_binding: Seq, n: int, rigour: int, sep: str = '\n',
                      V: bool = False, fasta: bool = False,
                      ind: int = 0) -> str:
    """Returns a list of <n> sets of primers, generated using the given
    components, with <rigour>. Seperates each set with <sep>. Will print
    additional information to console iff <V>. Will format data in fasta format
    iff fasta. Will shift the index of the printed primer by <ind> * <n> bases."""

    hg = HeteroGen(presenter=ConsolePresenter())

    set_list = []

    incomplete_forward_primer = MBPrimerBuilder(adapter_seq=for_adapter,
                                                binding_seq=for_binding)
    incomplete_reverse_primer = MBPrimerBuilder(adapter_seq=rev_adapter,
                                                binding_seq=rev_binding)

    # Choose n random spacers from the ten best.
    for_spacers = random.choices(
        hg.get_all_spacer_combos(incomplete_forward_primer._binding_seq)[
        0:n + 2],
        k=n)
    rev_spacers = random.choices(
        hg.get_all_spacer_combos(incomplete_reverse_primer._binding_seq)[
        0:n + 2],
        k=n)

    for set_num in range(n):

        # Select random spacer from the ten best
        for_spacer = for_spacers[set_num]
        rev_spacer = rev_spacers[set_num]

        hg.set_pairwise()
        hg.set_rigour(rigour)

        number_to_return = 1

        if V:
            print(for_spacer)
            print(rev_spacer)
            print("Generating primers. This may take some time...")

        # Generate the sets.
        t0 = time()
        psets = hg.get_hetero_seqs(incomplete_forward_primer,
                                   incomplete_reverse_primer,
                                   for_spacer, rev_spacer, number_to_return)
        runtime = time() - t0

        if V:
            print("Completed set #{num:d} in {time:.2f} seconds."
                  .format(time=runtime, num=set_num + n * ind, ))
        set_list.extend(psets)

    primer_str = ''
    for set in enumerate(set_list):
        if fasta:
            primer_str += set[1].get_fasta_seqs(set[0] + ind * n)
        else:
            primer_str += set[1].get_plain_seqs()
        primer_str += sep

    return (primer_str)


def get_primer_sets(gen_set: GenSet, n: int, rigour: int = RIGOUR,
                    filepath: Path = OUT_FILE,
                    sep: str = '\n', V: bool = VERBOSE, fasta: bool = DO_FASTA,
                    multi_file=MULTI_FILE,
                    filename: str = "Primer Set") -> None:
    """Runs gen_n_primers for each of the given binding sequences with the given
     parameters. Stores the output in text files in <filepath>. Will place
     all sets into a single file iff multi_file, othwise will create output file
     for each pair of forward and reverse primers."""
    if not gen_set.useable():
        raise ValueError("Gen set is unusable. Make sure it is properly "
                         "instantiated.")

    for_adapters, rev_adapters, for_bindings, rev_bindings = gen_set.unpack()
    out_str = ''
    for i in range(min(len(for_bindings), len(rev_bindings))):
        for_binding = for_bindings[i]
        rev_binding = rev_bindings[i]
        for_adapter = rev_adapters[i]
        rev_adapter = for_adapters[i]

        if not multi_file:
            ind = i

        # Get some primer sets as a string.
        s = get_n_primer_sets(for_adapter, rev_adapter, for_binding,
                              rev_binding, n, rigour, sep, V, fasta, ind=i)
        if fasta:
            txt_path = filepath / "{filename}@set{num:d}.fst".format(
                filename=filename,
                num=i + 1)
        else:
            txt_path = filepath / "{filename}@set{num:d}.txt".format(
                filename=filename,
                num=i + 1)
        # Store string in file.
        if multi_file and GEN_OUTPUT:
            with open(txt_path, "w") as my_file:
                my_file.write(s)
        else:
            out_str += s

    if not multi_file and GEN_OUTPUT:
        if fasta:
            txt_path = filepath / (filename + ".fst")
        else:
            txt_path = filepath / (filename + ".txt")
        if V:
            print(filename + " complete")
        with open(txt_path, "w") as my_file:
            my_file.write(out_str)


if __name__ == '__main__':
    main()
