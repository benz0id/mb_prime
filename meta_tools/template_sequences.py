from Bio.Seq import Seq
from typing import List

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
                    len(self.forward_adapters) == len(self.reverse_adapters):
                return True
            else:
                return False
        except AttributeError:
            return False

    def has_binding_num(self, ind: int) -> bool:
        return min(len(self.forward_binding), len(self.reverse_binding)) < ind


STANDARD_SET = GenSet()

STANDARD_SET.for_adapters = [Seq('ACACTCTTTCCCTACACGACGCTCTTCCGATCT')]
rev_adapters = [Seq('GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT')]

STANDARD_SET.for_bindings = [Seq('GTCAACCCAGCAGCTTACGCTGCCCTGGGTGCCTACA')]
STANDARD_SET.rev_bindings = [Seq('TGTGGCTTGGTGGATCTTCACACATCAGGGCTCTGAT')]


RANDOM_RHO_SET = GenSet()

RANDOM_RHO_SET.for_bindings = [Seq('GTCAACCCAGCAGCTTACGCTGCCCTGGGTGCCTACA'),
                Seq('TCACCACAACGATGTACACCTCCATGCAC'),
                Seq('CTTGGCTGCAATCTGGAAGGATTCTTTGC'),
                Seq('TCTCTGGTCACTGGTCGTTCTGGCTATTG'),
                Seq('AGCCCATCAGCAACTTCCGCTTCGGAGAG'),
                Seq('AACTACATCCTGCTGAATCTCGCGGTGGCCGACC')]
RANDOM_RHO_SET.rev_bindings = [Seq('TGTGGCTTGGTGGATCTTCACACATCAGGGCTCTGAT').reverse_complement(),
                Seq('TTTGCCAAGAGTTCCTCCATCTACAACCC').reverse_complement(),
                Seq('AGACCACCCAGAGGGCTGAGAGGGAAG').reverse_complement(),
                Seq('ATCATCATGGTCTTCGCCTTCCTGATG').reverse_complement(),
                Seq('CCGAGACCACCCAGAGGGCTGAGAGGGAA').reverse_complement(),
                Seq('AATTGTGTTCTTCTGCTACGGCCGTCTGCTGTGCGCC').reverse_complement()]


SMALL_BINDING_SET = GenSet()


SMALL_BINDING_SET.for_adapters = Seq('ACACT')
SMALL_BINDING_SET.rev_adapters = Seq('GTGAC')

SMALL_BINDING_SET.for_bindings = [Seq('GTCAA'),
                  Seq('TCACC'),
                  Seq('CTTGG'),
                  Seq('TCTCT'),
                  Seq('AGCCC'),
                  Seq('AACTA')]
SMALL_BINDING_SET.rev_bindings = [Seq('TGTGG').reverse_complement(),
                  Seq('TTTGC').reverse_complement(),
                  Seq('AGACC').reverse_complement(),
                  Seq('ATCAT').reverse_complement(),
                  Seq('CCGAG').reverse_complement(),
                  Seq('AATTG').reverse_complement()]
