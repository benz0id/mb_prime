from meta_tools.parse_tf_out import Homodimer, Heterodimer, parse_tf_output
from typing import List, Tuple, Dict
from pathlib import Path
import os
# Parses some input file from a TF run.



# Pre-set Parameter types parameters must be alpha and of length 3.
RIGOUR = 'rig'
SETNUM = 'set'
SETSIZE = 'sze'
REPLIC = 'rep'

FT = '.txt'
DELIM = '@'

class Parameter:
    """ A class designed to specify a parameter of a set of primers to be
    generated, or their computed properties.

    === Private Attributes ===
    param:
            The type of this parameter. One of a limited set of parameters
            types.

    val:
            The value of this parameter."""

    _param: str
    _val: int

    def __init__(self) -> None:
        """Initialises an empty Parameter"""
        return

    def set_param(self, param: str) -> None:
        self._param = param

    def set_val(self, val: int) -> None:
        self._val = val

    def get_param(self) -> str:
        return self._param

    def get_val(self) -> int:
        return self._val




def parse_parameters(filename: str, delim: str = DELIM,
                     file_type: str = FT) -> List[Parameter]:
    """Returns all parameters specified in a <filename>. Assumes that parameters
    are seperated by <delim>."""
    filename = filename.strip()
    received_filetype = filename[-len(file_type):]

    if file_type != received_filetype:
        raise ValueError("File must have filetype matching the one given.")

    params = []
    filename = filename[:-len(file_type)]
    raw_params = filename.split(delim)
    for raw_param in raw_params:
        param = Parameter()

        # Extract and check param type.
        param_type = raw_param[0:3]
        if len(param_type) > 3 or not param_type.isalpha():
            raise ValueError("Received invalid input parameter: " + raw_param)
        param.set_param(param_type)

        # Extract and check param value.
        param_val = raw_param[3:]
        try:
            param.set_val(int(param_val))
        except ValueError:
            raise ValueError("Received invalid input parameter: " + raw_param)
        params.append(param)

    return params


class Result:
    """The results from a single thermofischer primer analysis.

    === Private Attributes ===
    _parameters: The parameters with which this run was conducted.

    _homodimers: The homodimers found in this run.

    _heterodimers: The heterodimers found in this run."""

    _parameters: List[Parameter]
    _homodimers: List[Homodimer]
    _heterodimers:List[Heterodimer]

    def __init__(self, parameters: List[Parameter], homodimers: List[Homodimer],
                 heterodimers: List[Heterodimer]) -> None:
        """Initialises a Result with the given values."""
        self._parameters = parameters
        self._homodimers = homodimers
        self._heterodimers = heterodimers

    def __str__(self) -> str:
        """Returns a string representation of this result."""
        out_str = ''
        out_str += "Homodimers:\n"
        for homodimer in self._homodimers:
            out_str += '\t' + str(homodimer) + '\n'

        out_str += "Heterodimers:\n"
        for heterodimer in self._heterodimers:
            out_str += '\t' + str(heterodimer) + '\n'

        return out_str

    def optimise_dimers(self) -> None:
        """Removes all but the most stable homodimers for each conformation."""
        to_remove = []
        # For each unique combinations of homodimers, see if they specify the
        # same template primer.
        for i in range(len(self._homodimers)):
            for j in range(i + 1, len(self._homodimers)):
                h1 = self._homodimers[i]
                h2 = self._homodimers[j]
                same_primer = h1.get_primer_ind() == h2.get_primer_ind() and \
                              h1.get_is_forward() == h2.get_is_forward()
                # These primers are the same, so remove the least stable.
                if same_primer:
                    if h1.get_num_comp() < h2.get_num_comp():
                        to_remove.append(i)
                    else:
                        to_remove.append(j)
        # Ensure uniqueness.
        to_remove = list(set(to_remove))
        for ind in to_remove[::-1]:
            self._homodimers.pop(ind)

        to_remove = []
        # Repeat for heterodimers.
        for i in range(len(self._heterodimers)):
            for j in range(i + 1, len(self._heterodimers)):
                h1 = self._heterodimers[i]
                h2 = self._heterodimers[j]

                same_primer = h1.get_forward_ind() == h2.get_forward_ind() and \
                    h1.get_reverse_ind() == h2.get_reverse_ind()
                # These primers are the same, so remove the least stable.
                if same_primer:
                    if h1.get_num_comp() < h2.get_num_comp():
                        to_remove.append(i)
                    else:
                        to_remove.append(j)
        # Ensure uniqueness.
        to_remove = list(set(to_remove))
        for ind in to_remove[::-1]:
            self._heterodimers.pop(ind)


    def prune_mismatch_heteros(self) -> None:
        """Removes heterodimers that don't have matching forward and reverse
         indices."""
        to_remove = []
        for i in enumerate(self._heterodimers):
            hd = i[1]
            ind = i[0]
            if hd.get_forward_ind() != hd.get_reverse_ind():
                to_remove.append(ind)

        to_remove.sort(reverse=True)
        for ind in to_remove:
            self._heterodimers.pop(ind)

    def prune(self) -> None:
        """Removes superfluous dimers from this set."""
        self.prune_mismatch_heteros()
        self.optimise_dimers()

    def get_homodimers(self) -> List[Homodimer]:
        return self._homodimers

    def get_heterodimers(self) -> List[Heterodimer]:
        return self._heterodimers

    def get_parameters(self) -> List[Parameter]:
        return self._parameters

    def get_homo_scores(self) -> List[int]:
        """Returns the scores for the homodimers found in this run."""
        homo_scores = []
        for hd in self._homodimers:
            homo_scores.append(hd.get_num_comp())
        return homo_scores

    def get_hetero_scores(self) -> List[int]:
        """Returns the scores for the heterodimers found in this run."""
        hetero_scores = []
        for hd in self._heterodimers:
            hetero_scores.append(hd.get_num_comp())
        return hetero_scores

    def get_param_val(self, param_name: str) -> int:
        for param in self._parameters:
            if param.get_param() == param_name:
                return param.get_val()
        raise ValueError("Param does not exist, or is not a parameter of this "
                         "result.")


def get_all_results(dirpath: str or Path) -> List[Result]:
    """Parses a set of formatted TF results in <filepath> and stores them in a
    list of Result objects, which are then returned."""
    # Convert to path if neccesary
    if type(dirpath) == str:
        dirpath = Path(dirpath)

    filenames = os.listdir(dirpath)
    results = []
    # Compute result for each file in directory
    for filename in filenames:
        file_path = dirpath / filename
        file_params = parse_parameters(filename)
        file_homodimers, file_heterodimers = parse_tf_output(file_path)
        results.append(Result(file_params, file_homodimers, file_heterodimers))

    return results


def create_empty_output_files(parameters: Dict[str, Tuple[int, int]],
                              filepath: Path = '', filetype: str = FT,
                              make_files: bool = True) -> None or List[str]:
    """Creates empty text files in <filepath> with names given by <parameters>.
    A file name will be produced for every combination of parameters. Will not
    makes files iff not <make_files>, instead will return string of filenames.

    Precondition:
        make_files == True => <filepath> is a valid filepath.

    parameters: <str> being the name of the parameter and <Tuple[int, int]>
        being the range of values that parameter can span [int, int]. First
        index must be less than or equal to second index.


    >>> import os
    >>> filepath = Path("./Desktop") # an empty file
    >>> parameters = {RIGOUR: (1, 3), SETNUM: (1, 3)}
    >>> create_empty_output_files(parameters, filepath)
    >>>os.listdir(filepath)
    """
    # Calculate number of outfiles.
    ranges = parameters.values()
    num_outfiles = 1
    for range_inds in ranges:
        num_outfiles = num_outfiles * ( range_inds[1] - range_inds[0] + 1)

    # Generate empty outfiles names.
    filenames = []
    for i in range(num_outfiles):
        filenames.append('')

    num_params = len(parameters.keys())
    delim = DELIM
    repeat_length = num_outfiles
    # Append parameters to strings.
    for param in enumerate(parameters.keys()):
        range_inds = parameters[param[1]]
        param_name = param[1]
        param_ind = param[0]
        # We don't need a delimieter on the last param.
        if param_ind == num_params - 1:
           delim = ''

        num_repeats = num_outfiles // repeat_length
        # Append to filenames in tree-like fashion.
        repeat_length = repeat_length // (range_inds[1] - range_inds[0] + 1)
        for i in range(num_repeats):
            repeat_start = i * repeat_length * \
                           (range_inds[1] - range_inds[0] + 1)
            next_repeat_ind = 0
            # Each allowable value of a parameter according to the range.
            for param_num in range(range_inds[0], range_inds[1] + 1):
                last_repeat_ind = next_repeat_ind
                next_repeat_ind += repeat_length
                # Repeat to construct tree.
                param_str = param_name + str(param_num) + delim
                for i in range(last_repeat_ind, next_repeat_ind):
                    filenames[repeat_start + i] += param_str

    # Append filetype
    for i in range(len(filenames)):
        filenames[i] += filetype

    # We needn't make the files.
    if not make_files:
        return filenames

    # Make files.
    for filename in filenames:
        path = filepath / filename
        with open(path, 'w') as new_file:
            new_file.write('\n')

    return


def make_and_get_results(params: Dict[str, Tuple[int, int]], filepath: Path) \
        -> List[Result]:
    """Promps the user to make empty result files, then waits for them to enter
    data, and returns results."""

    inpt = input("Enter 'Y' if you'd like to create empty result files in " +
                 str(filepath) + '\n' + "Enter 'S' to skip this step.\n>")
    if inpt == 'S':
        pass
    elif inpt.strip() != 'Y':
        quit()

    create_empty_output_files(params, filepath)

    while input != 'R':
        inpt = input("Enter 'R' if ready to parse data in the above path or 'Q'"
                     " to quit.")
        if inpt == 'Q':
            quit()

    return get_all_results(filepath)




