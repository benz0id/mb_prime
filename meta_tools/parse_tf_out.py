from pathlib import Path
from typing import Tuple, List


def SeqMisalignedError(Exception):
    """Thrown when the parser becomes misaligned with the sequence."""
    pass


class Homodimer(object):
    """An object used to store the homodimer results of a TF analysis.

    === Private Attributes ===:
        _primer_ind: The index of the input primer that produced the homodimer.

        _is_forward: The direction of the primer

        _reverse_offset: By reversing this primer and shifting <reverse_offset>
            bases to the right, this particular homodimer is produced.

        _num_comp: The number of complementary bases in this homodimer.
    """

    _primer_ind: int
    _is_forward: bool
    _reverse_offset: int
    _num_comp: int

    def __init__(self) -> None:
        """Initialises an empty Homodimer object."""
        pass

    def __str__(self):
        """Returns a string representation of this heterodimer"""
        return "Homodimer {dir}{num: d}:Comp: {comp: d}, Offset: {off: d}".\
            format(
            dir = ['R', 'F'][self._is_forward], num = self._primer_ind,
            comp = self._num_comp, off = self._reverse_offset)

        # Getters and setters.
    def set_primer_ind(self, ind:int) -> None:
        self._primer_ind = ind
        return

    def set_is_forward(self, is_forward: bool) -> None:
        self._is_forward = is_forward
        return

    def set_reverse_offset(self, offset: int) -> None:
        self._reverse_offset = offset
        return

    def set_num_comp(self, num_comp: int) -> None:
        self._num_comp = num_comp
        return

    def get_primer_ind(self) -> int:
        return self._primer_ind

    def get_is_forward(self) -> bool:
        return self._is_forward

    def get_reverse_offset(self) -> int:
        return self._reverse_offset

    def get_num_comp(self) -> int:
        return self._num_comp

    def get_copy(self):
        """Returns a deep copy of this object"""
        hd = Homodimer()
        hd._num_comp = self._num_comp
        hd._is_forward = self._is_forward
        hd._primer_ind = self._primer_ind
        hd._reverse_offset = self._reverse_offset
        return hd


class Heterodimer:
    """An object used to store the heterodimer results of a TF analysis.

    === Private Attributes ===:
        _forward_ind: The index of the input forward primer.
        _reverse_ind : The index of the input reverse primer.

        _reverse_offset: By reversing the reverse primer (so that it is 3'-5'_
            and shifting <reverse_offset> bases to the right, this particular
            heterodimer is produced.

        _num_comp: The number of complementary bases in this homodimer.
    """
    _forward_ind: int
    _reverse_ind: int
    _reverse_offset: int
    _num_comp: int

    def __init__(self) -> None:
        """Initialises an empty Heterodimer object."""
        pass

    def __str__(self):
        """Returns a string representation of this heterodimer"""
        return "Heterodimer F{F}-R{R}. Comp: {comp: d}, Offset: {off: d}".format(
            F = self._forward_ind, R = self._reverse_ind, comp = self._num_comp,
            off = self._reverse_offset)

    # Getters and setters.
    def set_forward_ind(self, ind:int) -> None:
        self._forward_ind = ind
        return

    def set_reverse_ind(self, ind:int) -> None:
        self._reverse_ind = ind
        return

    def set_reverse_offset(self, offset: int) -> None:
        self._reverse_offset = offset
        return

    def set_num_comp(self, num_comp: int) -> None:
        self._num_comp = num_comp
        return

    def get_forward_ind(self) -> int:
        return self._forward_ind

    def get_reverse_ind(self) -> int:
        return self._reverse_ind

    def get_reverse_offset(self) -> int:
        return self._reverse_offset

    def get_num_comp(self) -> int:
        return self._num_comp

    def get_copy(self):
        """Returns a deep copy of this object"""
        hd = Heterodimer()
        hd._num_comp = self._num_comp
        hd._reverse_ind = self._reverse_ind
        hd._forward_ind = self._forward_ind
        hd._reverse_offset = self._reverse_offset
        return hd

def parse_tf_output(filepath: Path or str) -> Tuple[List[Homodimer],
                                                    List[Heterodimer]]:

    # Open file containing output from ThermoFisher multiple primer analyser.
    output_path = Path(Path(filepath))
    with open(output_path, 'r') as output_file:
        output_lines = output_file.readlines()

    # Will contain parsed homo/heterodimers
    homodimers = []
    heterodimers = []

    i = 0

    def cur_line() -> str:
        # Returns the current line according to i with newline removed.
        return output_lines[i][0:-1]

    mode = "und" # undefined
    while i < len(output_lines):

        # We've entered the self dimer region of the file.
        if "Self-Dimers:" in cur_line():
            mode = 'sd' # self dimer

        if "Cross Primer Dimers:" in cur_line():
            mode = "hd" # heterodimers

        # Parse homodimers.
        while mode == 'sd' and i < len(output_lines):

            if 'dimers for' in cur_line(): # we're on a line specifying some homodimers for a primer.

                # Collect attributes for this set of dimers.
                num_homodimers = int(cur_line()[0])
                primer_ind = int(cur_line()[-1])
                is_forward = cur_line()[-2] == 'F'
                i += 1
                template_dimer = Homodimer()
                template_dimer.set_is_forward(is_forward)
                template_dimer.set_primer_ind(primer_ind)

                # Parse each homodimer below.
                for _ in range(num_homodimers):
                    # We know that there's 4 lines in each homodimer.
                    homodimer_lines = output_lines[i: i + 4]
                    i += 4

                    # Quantify offset and complementarity.
                    comp = homodimer_lines[1].count('|')
                    offset = - homodimer_lines[0].count(' ') + \
                             homodimer_lines[2].count(' ')

                    # Catch misalignment
                    if homodimer_lines[3] != '\n':
                        raise SeqMisalignedError(
                            "Sequence misaligned at line {ind: d}".format(ind = i))

                    template_dimer.set_num_comp(comp)
                    template_dimer.set_reverse_offset(offset)
                    homodimers.append(template_dimer.get_copy())
            else:
                i += 1

            if "Cross Primer Dimers:" in cur_line():
                mode = "hd" # heterodimers


        while mode == 'hd' and i < len(output_lines):

            # At the start of a block specifying a
            if 'with' in cur_line():

                p1, p2 = cur_line()[:2], cur_line()[-2:]
                # Do we have two primers of the same kind? Record data if not.
                record = not p1[0] == p2[0]
                # Invert the reverse offset if necessary.
                invert = p1[0] == 'R'
                # Skip line showing leading primer. p1 is leading primer.
                i += 2

                # Find inds of forward and reverse primers.
                if p1[0] == 'F':
                    f_ind = int(p1[1])
                    r_ind = int(p2[1])
                else:
                    r_ind = int(p1[1])
                    f_ind = int(p2[1])

                # Format template dimer
                template_dimer = Heterodimer()
                template_dimer.set_forward_ind(f_ind)
                template_dimer.set_reverse_ind(r_ind)

                # Only parse  the homodimers for this pair of primers. Stop if we've
                # reached the end of the file or we're in the next set of primers.
                while i + 3 < len(output_lines) and (i + 4 > len(output_lines) - 1 or
                                                      "with" in output_lines[i + 4]):
                    # Read lines containing heterodimer data.
                    heterodimer_lines = output_lines[i:i+4]
                    i += 4
                    if invert:
                        # Calculate size of region beyond the 3' end of the rev
                        # primer.
                        reverse_offset = len(heterodimer_lines[2]) - \
                                         len(heterodimer_lines[0])
                    else:
                        reverse_offset = heterodimer_lines[2].count(' ')

                    comp = heterodimer_lines[1].count('|')
                    # Catch misalignment
                    if heterodimer_lines[3] != '\n':
                        raise SeqMisalignedError(
                            "Sequence misaligned at line {ind: d}".format(ind = i))

                    # Add dimer to list of dimers.
                    template_dimer.set_reverse_offset(reverse_offset)
                    template_dimer.set_num_comp(comp)
                    if record:
                        heterodimers.append(template_dimer.get_copy())
            else:
                i += 1
        i += 1
    return homodimers, heterodimers

