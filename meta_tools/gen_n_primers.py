from Bio.Seq import Seq

from hetero_spacer_generator.sequence_tools import is_degen
from presenters import ConsolePresenter
import demo_tools
from time import time
from hetero_spacer_generator.primer_tools import MBPrimerBuilder
from typing import Dict, List, Union, Tuple
from hetero_spacer_generator.defaults import ABSOLUTE_MAX_NUM_HETERO, \
    ABSOLUTE_MAX_SPACER_LENGTH
from time import sleep
import random
from hetero_spacer_generator.spacer_generator.hetero_spacer_generator import HeteroGen
from pathlib import Path


# Script to generate some number of primers for a given sequence

for_adapter = Seq('ACACTCTTTCCCTACACGACGCTCTTCCGATCT')
rev_adapter = Seq('GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT')

for_binding = Seq('GTCAACCCAGCAGCTTACGCTGCCCTGGGTGCCTACA')
rev_binding = Seq('TGTGGCTTGGTGGATCTTCACACATCAGGGCTCTGAT')

for_bindings = [Seq('GTCAACCCAGCAGCTTACGCTGCCCTGGGTGCCTACA'),
                Seq('TCACCACAACGATGTACACCTCCATGCAC'),
                Seq('CTTGGCTGCAATCTGGAAGGATTCTTTGC'),
                Seq('TCTCTGGTCACTGGTCGTTCTGGCTATTG'),
                Seq('AGCCCATCAGCAACTTCCGCTTCGGAGAG'),
                Seq('AACTACATCCTGCTGAATCTCGCGGTGGCCGACC')]
rev_bindings = [Seq('TGTGGCTTGGTGGATCTTCACACATCAGGGCTCTGAT').reverse_complement(),
                Seq('TTTGCCAAGAGTTCCTCCATCTACAACCC').reverse_complement(),
                Seq('AGACCACCCAGAGGGCTGAGAGGGAAG').reverse_complement(),
                Seq('ATCATCATGGTCTTCGCCTTCCTGATG').reverse_complement(),
                Seq('CCGAGACCACCCAGAGGGCTGAGAGGGAA').reverse_complement(),
                Seq('AATTGTGTTCTTCTGCTACGGCCGTCTGCTGTGCGCC').reverse_complement()]

for_adapter_s = Seq('ACACT')
rev_adapter_s = Seq('GTGAC')

for_binding_s = Seq('GTCAACCCAGCAGCTTACGCTGCCCTGGGTGCCTACA')
rev_binding_s = Seq('TGTGGCTTGGTGGATCTTCACACATCAGGGCTCTGAT')

for_bindings_s = [Seq('GTCAA'),
                Seq('TCACC'),
                Seq('CTTGG'),
                Seq('TCTCT'),
                Seq('AGCCC'),
                Seq('AACTA')]
rev_bindings_s = [Seq('TGTGG').reverse_complement(),
                Seq('TTTGC').reverse_complement(),
                Seq('AGACC').reverse_complement(),
                Seq('ATCAT').reverse_complement(),
                Seq('CCGAG').reverse_complement(),
                Seq('AATTG').reverse_complement()]

class GenSet:
    """Stores primers to be used in testing of primer generation algorithms.

    All attributes are public.
    """
    forward_binding: List[Seq]
    reverse_binding: List[Seq]
    forward_adapters: List[Seq]
    reverse_adapters: List[Seq]
    def __init__(self) -> None:
        """Initialises and empty GenSet"""
        return

    def useable(self) -> bool:
        """Returns whether this class has all attributes instantiated with
        proper values."""
        try:
            if len(self.forward_binding) == len(self.reverse_binding) and \
                len(for_adapter) == len(rev_adapter):
                return True
            else:
                return False
        except AttributeError:
            return False

    def has_binding_num(self, ind: int) -> bool:
        return min(len(self.forward_binding), len(self.reverse_binding)) < ind


RIGOUR = 1
NUM_SETS = 3
OUT_FILE = Path('C:\\Users\\bfern\\OneDrive - University of Toronto\\Desktop\\Short Primers')
VERBOSE = True
DO_FASTA = True


def main():
    for rigour in range(0, 3):
        filename = "Rigour " + str(rigour)
        print(filename)
        get_primer_sets(for_adapter, rev_adapter, for_bindings, rev_bindings, NUM_SETS,
                    rigour=rigour, filepath=OUT_FILE, V=VERBOSE, fasta=DO_FASTA,
                        multi_file=False, filename=filename)


def get_n_primer_sets(for_adapter: Seq, rev_adapter: Seq, for_binding: Seq,
                  rev_binding : Seq, n: int, rigour: int, sep: str = '\n',
                  V: bool = False, fasta: bool = False, ind: int  = 0) -> str:
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
    for_spacers = random.choices(hg.get_all_spacer_combos(incomplete_forward_primer._binding_seq)[0:n + 2],
                                k = n)
    rev_spacers = random.choices(hg.get_all_spacer_combos(incomplete_reverse_primer._binding_seq)[0:n + 2],
                                k = n)

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
        psets = hg.get_hetero_seqs(incomplete_forward_primer, incomplete_reverse_primer,
                                   for_spacer, rev_spacer, number_to_return)
        runtime = time() - t0

        if V:
            print("Completed set #{num:d} in {time:.2f} seconds."
                  .format(time=runtime, num = set_num + n * ind,))
        set_list.extend(psets)


    primer_str = ''
    for set in enumerate(set_list):
        if fasta:
            primer_str += set[1].get_fasta_seqs(set[0] + ind * n)
        else:
            primer_str += set[1].get_plain_seqs()
        primer_str += sep

    return(primer_str)

def get_primer_sets(for_adapter: Seq, rev_adapter: Seq, for_bindings: List[Seq],
                    rev_bindings: List[Seq], n: int, rigour: int = RIGOUR, filepath: Path = OUT_FILE,
                    sep: str = '\n', V: bool = False, fasta: bool = False,
                    multi_file = True, filename: str = "Primer Set") -> None:
    """Runs gen_n_primers for each of the given binding sequences with the given
     parameters. Stores the output in text files in <filepath>. Will place
     all sets into a single file iff multi_file, othwise will create output file
     for each pair of forward and reverse primers."""


    out_str = ''
    for i in range(min(len(for_bindings), len(rev_bindings))):
        for_binding = for_bindings[i]
        rev_binding = rev_bindings[i]

        if not multi_file:
            ind = i

        # Get some primer sets as a string.
        s = get_n_primer_sets(for_adapter, rev_adapter, for_binding,
        rev_binding, n, rigour, sep, V, fasta, ind=i)
        if fasta:
            txt_path = filepath / "{filename}@set{num:d}.fst".format(filename = filename,
                num = i + 1)
        else:
            txt_path = filepath / "{filename}@set{num:d}.txt".format(filename = filename,
            num = i + 1)
        # Store string in file.
        if multi_file:
            with open(txt_path, "w") as my_file:
                my_file.write(s)
        else:
            out_str += s

    if not multi_file:
        if fasta:
            txt_path = filepath / (filename + ".fst")
        else:
            txt_path = filepath / (filename + ".txt")
        if V:
            print(filename + " complete")
        with open(txt_path, "w") as my_file:
            my_file.write(out_str)


def gen_n_by_param(parameters: Dict[str, Tuple[int, int]],
                   template_seqs: GenSet) -> None:
    """Generates """


if __name__ == '__main__':
    main()
