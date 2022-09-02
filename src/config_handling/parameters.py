# functions for extracting and storing config files.
from typing import Any, Dict, List, Tuple
import logging
import abc
import os

from src.config_handling import input_validator as iv, command_line_tools as cli, \
    formatting as fmt
from pathlib import Path
from test_files.fixtures_and_helpers import ALIGNMENTS_PATH

from src.seq_alignment_analyser.align import MSA

log = logging.getLogger('root')


class Parameter(abc.ABC):
    """A parameter to be stored in the config file by the user. Used to query
    the user and generate a string that can be stored in a python file.
    === Public Attributes ===
    name: The name of this parameter."""

    name: str

    def __init__(self, name: str) -> None:
        """Initializes the class using the given values."""
        self.name = name

    @abc.abstractmethod
    def query_user(self) -> None:
        """Queries the user to provide the value of this parameter."""
        pass

    @abc.abstractmethod
    def get_python_variable_string(self) -> str:
        """Returns a string representation of this variable, ready to be pasted
        into a file as a valid variable."""
        pass


class SimpleParameter(Parameter):
    """
    A parameter to be stored in the config file by the user. Used to query the
    user and generate a string that can be stored in a python file.

    === Private Attributes ===
    _validations: Function used to validate the input when received from user.
    _err_msg: To be printed when receiving invalid input.

    === Public Attributes ===
    data: The value of this parameter.
    """

    _validations: List[iv.Validation]
    _msg: str
    data: Any

    def __init__(self, name: str, msg: str,
                 validations: List[iv.Validation]) -> None:
        """Initialises this parameter using the given values."""
        super().__init__(name)
        self._validations = validations
        self._msg = msg
        self.query_user()


class StrParam(SimpleParameter):
    """Responsible for storing and fetching string parameters.

    === Public Variables ===
    data: The string value of this parameter.
    """
    data: str

    def __init__(self, name: str, msg: str, validations: List[
        iv.Validation]
                 ) -> None:
        """Initialises the variable to the given values."""
        super().__init__(name, msg, validations)

    def get_python_variable_string(self) -> str:
        """Returns the formatted string value."""
        return self.name + ' = ' + repr(self.data)

    def query_user(self) -> None:
        """Fetch an integer value from the user."""
        s = cli.prompt(self._msg, self._validations)
        self.data = s


class IntParam(SimpleParameter):
    """Responsible for storing and fetching integer parameters.

    === Public Variables ===
    data: The integer value of this parameter.
    """
    data: int

    def __init__(self, name: str, msg: str, validations: List[
        iv.Validation]) -> None:
        """Initialises the variable to the given values."""
        valids = validations[:]
        valids.insert(0, iv.VALID_INT)
        super().__init__(name, msg, valids)

    def get_python_variable_string(self) -> str:
        """Returns the formatted integer value."""
        return self.name + ' = ' + str(self.data)

    def query_user(self) -> None:
        """Fetch an integer value from the user."""
        s = cli.prompt(self._msg, self._validations)
        self.data = int(s)


class FloatParam(SimpleParameter):
    """Responsible for storing and fetching integer parameters.

    === Public Variables ===
    data: The integer value of this parameter.
    """
    data: float

    def __init__(self, name: str, msg: str, validations: List[
        iv.Validation]) -> None:
        """Initialises the variable to the given values."""
        valids = validations[:]
        valids.insert(0, iv.VALID_FLOAT)
        super().__init__(name, msg, valids)

    def get_python_variable_string(self) -> str:
        """Returns the formatted integer value."""
        return self.name + ' = ' + str(self.data)

    def query_user(self) -> None:
        """Fetch an integer value from the user."""
        s = cli.prompt(self._msg, self._validations)
        self.data = float(s)



class TimeParam(SimpleParameter):
    """Responsible for storing and fetching time parameters.

    === Public Variables ===
    data: The integer value of this parameter.
    """
    data: fmt.TimeSpec
    hours: int
    minutes: int
    seconds: int

    def __init__(self, name: str, msg: str, validations: List[
        iv.Validation]) -> None:
        """Initialises the variable to the given values."""
        valids = validations[:]
        valids.insert(0, iv.VALID_TIME)
        super().__init__(name, msg, valids)

    def get_python_variable_string(self) -> str:
        """Returns the formatted integer value."""
        return self.name + ' = ' + str(self.data)

    def query_user(self) -> None:
        """Fetch an integer value from the user."""
        s = cli.prompt(self._msg + fmt.TIME, self._validations)
        self.data = fmt.TimeSpec(*fmt.to_time_spec(s))


class SeqPairParam(SimpleParameter):
    """Responsible for storing and fetching integer parameters.

    === Public Variables ===
    data: The value of this parameter.
    start: Start point of this range.
    stop: Stop point of this range.
    """
    data: fmt.SeqPairParam
    forward: str
    reverse: str

    def __init__(self, name: str, msg: str, validations: List[
        iv.Validation]) -> None:
        """Initialises the variable to the given values."""
        valids = validations[:]
        valids.insert(0, iv.VALID_DNA)
        super().__init__(name, msg, valids)

    def __str__(self) -> str:
        return ''.join(['SeqPairParam(', repr(self.forward), ', ',
                        repr(self.reverse), ')'])

    def get_python_variable_string(self) -> str:
        """Returns the formatted integer value."""
        return ''.join([
            self.name, ' = SeqPairParam', repr(self.data),
        ])

    def query_user(self) -> None:
        """Fetch an integer value from the user."""
        print(self._msg)
        self.forward = cli.prompt('Enter Forward Sequence 5\'-3\':',
                                  self._validations)
        self.reverse = cli.prompt('Enter Reverse Sequence 5\'-3\':',
                                  self._validations)


class RangeParam(SimpleParameter):
    """Responsible for storing and fetching integer parameters.

    === Public Variables ===
    data: The value of this parameter.
    start: Start point of this range.
    stop: Stop point of this range.
    """
    data: fmt.InclRange[int, int]
    start: int
    stop: int

    def __init__(self, name: str, msg: str, validations: List[
        iv.Validation]
                 , data: fmt.InclRange = None) -> None:
        """Initialises the variable to the given values."""
        valids = validations[:]
        valids.insert(0, iv.VALID_RANGE)
        super().__init__(name, msg, valids)
        if data:
            self.data = data

    def get_python_variable_string(self) -> str:
        """Returns the formatted integer value."""
        return ''.join([
            self.name, ' = ', str(self.data)
        ])

    def query_user(self) -> None:
        """Fetch an integer value from the user."""
        s = cli.prompt(self._msg + fmt.RANGE, self._validations)
        self.data = fmt.InclRange(*fmt.to_range(s))
        self.start, self.stop = self.data


class PathParam(SimpleParameter):
    """Responsible for storing and fetching integer parameters.

    === Public Variables ===
    data: The value of this parameter.
    """
    data: Path
    default: Path

    def __init__(self, name: str, msg: str, validations: List[
        iv.Validation], data: Path = None, default: Path = None) -> None:
        """Initialises the variable to the given values."""
        valids = validations[:]
        valids.insert(0, iv.VALID_PATH)
        if data:
            self.data = data
        if default:
            self.default = default
        else:
            self.default = Path(os.path.dirname(__file__)).parent
        super().__init__(name, msg, valids)

    def get_python_variable_string(self) -> str:
        """Returns the formatted integer value."""
        return ''.join([
            self.name, ' = Path(\'', str(self.data), '\')'
        ])

    def query_user(self) -> None:
        """Fetch an integer value from the user."""
        s = cli.prompt(self._msg, self._validations)
        if 'DIR' in s:
            self.data = self.default
        else:
            self.data = Path(s)


class TargetRegionParam(Parameter):
    """Responsible for storing and fetching integer parameters.

    === Private Variables ===

    _other_targets: To keep track of sites on this alignment already being
        targeted by other primer sets. Used to check for overlap between targets
        and/or primers. Also used to check that this target has a unique name.

    _sister_targets: Targets that lie on the same alignment as this one.

    _target_msa: The MSA this targets.

    _start_end_distance: The minimum distance between the primer and the start/
        end of alignment.

    === Public Variables ===
    data: The integer value of this parameter.

    filename: The name of the alignment these targets are targeting.

    min_sep_distance: The minimum distance between new primer sites and any used
        primer sites. To Account for primer binding sites.


    target_number: To help the user keep track of how many other sets of targets
        lie on this alignment.
    """
    data: fmt.TargetRegionInfo
    filename: str
    name: str
    target_sites: Tuple[int, ...]

    # Attributes for validating user input.
    min_sep_distance: int
    target_number: int
    _start_end_distance: int
    _other_targets: Dict[MSA, List[fmt.TargetRegionInfo]]
    _sister_targets: List[fmt.TargetRegionInfo]
    _target_msa: MSA

    def __init__(self, filename: str, max_target_len: int,
                 min_sep_distance: int, target_number: int,
                 start_end_distance: int,
                 other_targets: Dict[MSA, List[fmt.TargetRegionInfo]]) -> None:
        """Initialises the variable to the given values.

        Note:
            """
        super().__init__('Unnamed target config.')
        self.target_number = target_number
        self.max_target_length = max_target_len
        self.min_sep_distance = min_sep_distance
        self.filename = filename
        self._other_targets = other_targets
        self._start_end_distance = start_end_distance

        found = False
        for msa in other_targets.keys():
            if msa.filename == filename:
                self._target_msa = msa
                found = True
                break
        if not found:
            raise ValueError('Target MSA does not exist.')
        self._sister_targets = other_targets[self._target_msa]

    def __str__(self) -> str:
        return str(self.data)

    def get_python_variable_string(self) -> str:
        """Returns the formatted integer value."""
        return str(self)

    def _get_name(self) -> None:
        """Prompts the user to name this target."""
        used_names = []
        for targets in self._other_targets.values():
            for target in targets:
                used_names.append(target.name)

        vali = iv.Validation(lambda s: s not in used_names,
                             'That name has already been used.')

        self.name = cli.prompt('Enter a name for this set of targets.', [vali])

    def _get_target_sites(self) -> None:
        """Gets the sites that should be targeted in this set."""

        is_valid_format = iv.Validation(iv.is_list_of_ints, 'Input format: ' +
                                        str(fmt.INT_LIST))

        existing_targets = []
        for sister_target in self._sister_targets:
            existing_targets.extend(sister_target.sites)

        max_site_p1 = len(self._target_msa) - self._start_end_distance
        in_alignment = lambda s: iv.all_in_range(s, self._start_end_distance,
                                                 max_site_p1)

        sites_in_al = iv.Validation(in_alignment, 'Targets do not lie a '
                                                  'reasonable distance from the start and end of the alignment.')

        not_in_existing = iv.Validation(lambda s: iv.all_not_in(s,
                                                                existing_targets),
                                        'Contains sites already being '
                                        'targeted.',
                                        is_warning=True)

        space_between_existing = lambda s: iv.distance_gt(s, existing_targets,
                                                          self.min_sep_distance)
        err_msg = 'These targets overlap with another set that you have input.'
        adequate_space = iv.Validation(space_between_existing, err_msg,
                                       is_warning=True)

        all_within_target = lambda s: iv.internal_distance_lt(s,
                                                              self.max_target_length - 1)
        all_within_target = iv.Validation(all_within_target,
                                          'The specified sites are too far '
                                          'apart to fit into a single read.')

        s = cli.prompt('Enter any number of nucleotide positions you would '
                       'like to include in this target. ' + fmt.INT_LIST,
                       [
                           is_valid_format, sites_in_al, not_in_existing,
                           adequate_space, all_within_target
                       ])

        self.target_sites = tuple(fmt.to_list_of_ints(s))

    def query_user(self) -> None:
        """Fetch value from the user."""
        cli.print_title('Configuration for target set #' +
                        str(self.target_number) + ' on alignment '
                        + self.filename)

        self._get_name()

        self._get_target_sites()

        self.data = fmt.TargetRegionInfo(self.name, self.filename,
                                     self.target_sites)

        self._sister_targets.append(self.data)


def format_as_pylist(params: List[Parameter], name: str) -> str:
    """Produces string that can be read as a python variable using the given
    params. This variable will be of the format List[TargetRegionInfo]."""

    param_strs = [str(param) for param in params]

    list_contents = '\t' + ',\n\t'.join(param_strs)

    formatted_list = '[\n' + list_contents + '\n\t]'

    complete_list = name + ' = ' + formatted_list

    return complete_list


def simple_config() -> None:
    filename = 'alex_test_align.fas'
    aln_path = ALIGNMENTS_PATH / filename
    aln = MSA(aln_path)
    max_target_len = 150
    min_sep_distance = 40
    min_primer_size = 20
    targets_dict = {aln: []}
    target_nmbr = 1
    another = True
    while another:
        new_param = TargetRegionParam(filename, max_target_len,
                                      min_sep_distance, target_nmbr,
                                      min_primer_size, targets_dict)
        new_param.query_user()
        target_nmbr += 1
        another = cli.yes_no_prompt('Make another target?')
    print(format_as_pylist(targets_dict[aln], 'variable_name'))


if __name__ == '__main__':
    simple_config()
