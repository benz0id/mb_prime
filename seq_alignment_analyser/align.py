import os
from typing import Any, Dict, List, Tuple, Union
from pathlib import Path
from hetero_spacer_generator.defaults import DEGEN_TO_POSSIBLE, BASES, MIN_COMP, \
    WINDOW_SIZE
from hetero_spacer_generator.sequence_tools import is_degen_base
from Bio import SeqIO, AlignIO
import matplotlib.pyplot as plt
import logging

log = logging.getLogger('root')
DEFAULT_WINDOW_SIZE = 5


# An adapter class designed to parse alignments, and extract properties from
# them relevant to metabarcoding.
def add(dic: Dict, key: Any, item: Any) -> None:
    """Adds the <item> to the <dic> with <key>. Adds a list at <key> if it does
    not already exist in <dic>."""
    if key in dic:
        dic[key].append(item)
    else:
        dic[key] = [item]


def get_complementarity(b1: str, b2: str) -> float:
    """Returns the complementarity between the two given bases.
    EX A-T is complementary 100% of the time, A-N is complementary 1/4 of the
    time, and N-N is complementary """
    if not is_degen_base(b1) and not is_degen_base(b2):
        return b1 == b2
    b1_pos = DEGEN_TO_POSSIBLE[b1]
    b2_pos = DEGEN_TO_POSSIBLE[b2]
    total_perms = len(b1_pos) * len(b2_pos)
    equal_perms = 0
    for pb1 in b1_pos:
        for pb2 in b2_pos:
            equal_perms += pb1 == pb2

    return equal_perms / total_perms


def sliding_window(arr: List[Union[int, float]],
                   window_size: int) -> List[float]:
    """Returns a sliding distraction of the given array with section sizes of
    <window_size>."""
    windowed = []
    for i in range(len(arr) - window_size):
        avg = sum(arr[i: i + window_size]) / window_size
        windowed.append(avg)

    return windowed


KNOWN_MSA_TYPES = ['clustal', 'emboss', 'fasta', 'fasta-m10', 'ig',
                   'maf', 'mauve', 'msf', 'nexus', 'phylip',
                   'phylip-sequential', 'phylip-relaxed', 'stockholm']


class MSA:
    """
    A multiple sequence alignment containing a region to be targeted for
    metabarcoding analyses.
    === Public Attributes ===

    filename: The name of this MSA's containing file.

    === Private Attributes ===

    _seqs:
        An array storing the alignment. [id][base_pos] => base
    _seq_names:
        An array storing the alignment. [id] => name
    _window_size:
        The size of the distraction when generating sliding distraction graphs.
    _consensus:
        The consensus sequence of the alignment.
    _conservation:
        Calculated as the % of nucleotides at that position that are equal to
        the consensus sequence. Spacers are not counted.
    _percent_spacer:
        What percent of the sequences have a spacer at this position.
    _percent_missed:
        The % of bases at this position that we would miss with primer as they
        contain spacers, or do not contain spacers (when % spacers is over 50%).
        Never greater than 50%, as we would design primers to skip this region
        if only a few sequences had indels at this location.
    """

    filename: str
    _seqs: List[str]
    _seq_names: List[str]
    _window_size: int
    _consensus: str
    _conservation: List[float]
    _percent_spacer: List[float]
    _percent_missed: List[float]

    def __init__(self, filepath: Path, window_size: int = 5,
                 filetype: str = 'fasta'):
        """Initialises this MSA using information sored in <filepath>."""
        self.filename = os.path.basename(filepath)
        self._window_size = window_size
        self._parse_MSA(filepath, filetype)
        self._parse_consensus_attributes()

    def __hash__(self) -> int:
        """Hash is dependent on filename."""
        return hash(self.filename)

    def _parse_MSA(self, filepath: Path, filetype: str) -> None:
        """Initialises the sequence attributes of this class using the data
        stored in <filepath>. Accepts several filetypes."""
        self._seqs = []
        self._seq_names = []
        file_type = str(filepath).split('.')[-1]
        if file_type in ['fas', 'fst', 'fasta']:
            self._parse_from_fasta(filepath)
        else:
            self._parse_from_other(filepath, filetype)

        log.info(''.join([str(filepath), ' successfully parsed as a ',
                          filetype, ' file.']))
        return

    def _parse_from_other(self, filepath: Path, filetype: str) -> None:
        """Parses a non-fasta alignment."""
        if filetype not in KNOWN_MSA_TYPES:
            log.critical('Unknown file type received.')
            raise ValueError(''.join([
                '"', filetype, '"', 'is an unknown MSA filetype.',
                ' Known alignments include:', '\n'.join(KNOWN_MSA_TYPES)
            ]))
        try:
            alignment = AlignIO.read(open(filepath), filetype)
        except ValueError:
            log.critical('Failed to parse ' + str(filepath) + ' as a ',
                         filetype, ' file.')
            raise (IOError('Failed to parse ' + str(filepath) + ' as a ',
                           filetype, ' file.'))

        for record in alignment:
            self._seqs.append(record.seq)
            self._seq_names.append(record.id)

    def add_primers_to_graph(self, f_region_start: int, f_region_end: int,
                             r_region_start: int, r_region_end: int,
                             height: int = 50,
                             primer_region_name: str = 'Primer Regions',
                             amplicon_region_name: str = 'Amplicon Region') \
            -> None:
        """Highlights the amplicon_regions and primer regions on the graph. Sets
        the amplicon_regions to appear at <height>."""
        primer_regions = []
        p_data = []
        # Just add a bunch of blank data at height
        for i in range(f_region_start, f_region_end + 1):
            primer_regions.append(i)
            p_data.append(height)
        for i in range(r_region_start, r_region_end + 1):
            primer_regions.append(i)
            p_data.append(height)

        amplicon_regions = []
        a_data = []
        for i in range(f_region_end + 1, r_region_start):
            amplicon_regions.append(i)
            a_data.append(height)

        plt.plot(primer_regions, p_data, c='r', label=primer_region_name)
        plt.plot(amplicon_regions, a_data, c='g', label=amplicon_region_name)

        plt.legend(loc="lower left")

    def gen_plot(self) -> None:
        """Generates a plot showing the key attributes of this MSA."""
        plt.figure(figsize=(int(self.__len__() / 10), 6), dpi=80)

        x = list(range(1, len(self) + 1 - self._window_size))
        plt.tick_params(axis='x', which='major', labelsize=8)

        plt.title('Alignment Properties')

        conservation = sliding_window(self._conservation, self._window_size)
        percent_missed = sliding_window(self._percent_missed, self._window_size)

        plt.plot(x, conservation, c='g', label='Conservation')
        plt.plot(x, percent_missed, c='y', label='Sequences Missed')

        plt.legend(loc="upper left")

        fontsize = 10
        plt.xlabel(
            'Base Position (Window to i + ' + str(self._window_size) + ')',
            fontsize=fontsize)
        plt.ylabel('Percent of Sequences', fontsize=fontsize)

    def show_plot(self) -> None:
        """Shows the currently stored plot."""
        plt.show()

    def save_plot(self, filename: Union[Path, str]) -> None:
        """Saves the path using the name and location specified in <filename>.
        """
        x = list(range(1, len(self) + 1 - self._window_size))
        plt.xticks(x[4::5])
        plt.savefig(filename)

    def _parse_from_fasta(self, filepath: Union[Path, str]):
        """Parses MSA attributes from <filepath>."""
        fasta_sequences = SeqIO.parse(open(filepath), 'fasta')
        for fasta in fasta_sequences:
            name, sequence = fasta.description, str(fasta.seq)
            self._seqs.append(sequence)
            self._seq_names.append(name)

    def _parse_consensus_attributes(self) -> None:
        """Parses consensus attributes such as the consensus sequence and
        conservation from currently stored <self._seqs>.

        Precondition:
            There must be more of any one of the four bases than degenerate
            bases at any point in the alignment.
        """
        self._consensus = ''
        self._conservation = []
        self._percent_spacer = []
        self._percent_missed = []

        aln_len = len(self._seqs[0])
        num_seqs = len(self._seqs)

        # Check all inputs sequences are the same length.
        for i, seq in enumerate(self._seqs):
            if len(seq) != aln_len:
                raise ValueError("All input sequences must be of the same "
                                 "length. Offending sequence is #" + str(i) +
                                 ', ' + self._seq_names[i] + '.')

        # Extract attributes at each position.
        for base_ind in range(aln_len):

            # Find consensus.
            base_counts = {'A': 0, 'T': 0, 'C': 0, 'G': 0, '-': 0, 'other': 0}
            for seq in self._seqs:
                if seq[base_ind] not in 'ATCG-':
                    base_counts['other'] += 1
                else:
                    base_counts[seq[base_ind]] += 1

            # Extract most common base
            cons_base = max(BASES, key=lambda b: base_counts[b])

            # Append attributes for this base to lists of attributes.
            self._consensus += cons_base
            self._conservation.append(base_counts[cons_base] / num_seqs * 100)
            ps = base_counts['-'] / num_seqs * 100
            self._percent_spacer.append(ps)
            self._percent_missed.append(min(ps, abs(ps - 100)))

    def get_seqs(self) -> Dict[str, str]:
        """Returns a dictionary containing the sequences in this alignment,
        [name: sequence]."""
        seqs = {}
        for i, seq in enumerate(self._seqs):
            seqs[self._seq_names[i]] = seq
        return seqs

    def get_conservation_arr(self) -> List[float]:
        """Returns the conservation at the given <base_ind>."""
        return self._conservation[:]

    def get_mean_conservation(self, start: int, stop: int) -> float:
        """Returns the mean % conservation across the given region -
        [start, stop)."""
        s = 0
        for cons in self._conservation[start: stop + 1]:
            s += cons
        return s / (stop - start + 1)

    def get_percent_spacer_arr(self) -> List[float]:
        """Returns the percent of the alignment that is spacer at the given
        <base_ind>."""
        return self._percent_spacer[:]

    def get_percent_missed_arr(self) -> List[float]:
        """Returns the percent of the alignment that would be made out of reach
        of sequencing by having a primer bind at the given <base_ind>."""
        return self._percent_missed[:]

    def get_consensus(self) -> str:
        """Returns the consensus sequence of this alignment."""
        return self._consensus

    def get_conservation(self, base_ind: int) -> float:
        """Returns the conservation at the given <base_ind>."""
        return self._conservation[base_ind]

    def get_percent_spacer(self, base_ind: int) -> float:
        """Returns the percent of the alignment that is spacer at the given
        <base_ind>."""
        return self._percent_spacer[base_ind]

    def get_percent_missed(self, base_ind: int) -> float:
        """Returns the percent of the alignment that would be made out of reach
        of sequencing by having a primer bind at the given <base_ind>."""
        return self._percent_missed[base_ind]

    def __len__(self) -> int:
        """Returns the number of nucleotides in this alignment."""
        return len(self._seqs[0])

    def get_num_seqs(self) -> int:
        """Returns the number of sequences in this alignment."""
        return len(self._seqs)

    def scan_region(self, start: int, stop: int) -> str:
        """Scans the given region for possible complications that may arise by
        the placement of a primer in that location."""

        if not 0 <= start < stop < len(self._seqs[0]):
            raise ValueError("Invalid start and stop. Ensure the specified "
                             "region lies within this alignment.")

        species_to_missed_p = {}
        # Find polymorphisms that may be lost.
        for base_ind in range(start, stop):
            for i, seq in enumerate(self._seqs):
                if is_degen_base(seq[base_ind]):
                    add(species_to_missed_p, self._seq_names[i], base_ind)

        # Find max potential for loss of sequences.
        max_lost_ind = max(list(range(start, stop)),
                           key=lambda a: self.get_percent_missed(a))
        max_loss = self.get_percent_missed(max_lost_ind)

        low_cons = []
        # Find sites with low conservation.
        for base_ind in range(start, stop):
            if self.get_conservation(base_ind) < MIN_COMP:
                low_cons.append(base_ind)

        # Convert the above warnings into a string.
        s = ''
        if len(species_to_missed_p.keys()) > 0:
            s += ''.join([
                str(len(species_to_missed_p.keys())),
                ' species have polymorphisms lying in the selected region, '
                'including:\n'
            ])
            l = []
            for key in species_to_missed_p.keys():
                l.extend(['    ', key, ', at residue(s) ',
                          str(species_to_missed_p[key])[1:-1], '\n'])
            s += ''.join(l)
            s += '\n'

        if max_loss > 0:
            s += 'Due to indels at position ' + str(max_lost_ind) + \
                 ', a minimum of ' + str(
                max_loss) + ' % of target sequences will ' \
                            'experience severe binding disruption.\n\n'

        if len(low_cons) > 0:
            s += 'The following site have conservation levels below the minimum' \
                 ' threshold (' + str(MIN_COMP) + '%):\n'
            for ind in low_cons:
                s += '    ' + str(ind) + ':' + \
                     str(self.get_conservation(ind)) + '%\n'

        if s:
            a = 'Reporting warnings found in the selected binding region...\n\n'
            a += s
            return a
        return ''
