from typing import Any, Callable, List, Type, Union
from pathlib import Path
import logging
import sys
import config_handling.formatting as fmt
import seq_alignment_analyser.align as align

ADD_FORMAT_STR = True

log = logging.getLogger('root')


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


class Validation:

    method: Callable[[str], bool]
    err_msg: str
    is_warning: bool

    def __init__(self, method: Callable[[str], bool], msg: str,
                 is_warning: bool = False) -> None:
        self.method = method
        self.err_msg = msg
        self.is_warning = is_warning

    def __str__(self) -> str:
        return self.err_msg

    def is_valid(self, s: str) -> bool:
        return self.method(s)


# === Input validation functions ===


def in_range_incl(x: Union[int, str], start: int, stop: int) -> bool:
    """Returns whether <x> is in the given range."""
    if isinstance(x, str):
        x = int(x)
    return start <= x <= stop

def valid_int(s: str) -> bool:
    """Returns whether <s> can be interpreted as an integer."""
    try:
        fmt.int_strip(s)
    except ValueError:
        return False
    return True


def valid_float(s: str) -> bool:
    """Returns whether <s> can be interpreted as a float."""
    try:
        f = float(s)
    except ValueError:
        return False
    return True


def distance_gt(group1: str, group2: List[int], min_distance: int) -> bool:
    """Given that group1 is a validly formatted list of integers, tests that
    the difference between any two values in <group1> and <group2> is no less
    than min_distance."""
    group1 = fmt.to_list_of_ints(group1)
    for val1 in group1:
        for val2 in group2:
            if abs(val1 - val2) < min_distance:
                return False
    return True


def all_in_range(s: str, start: int, stop: int) -> bool:
    """Returns whether all of the values in <s> lie within start and stop."""
    ints = fmt.to_list_of_ints(s)
    for val in ints:
        if not start < val < stop:
            return False
    return True


def all_not_in(s: str, group: List[int]) -> bool:
    """Given <s> is a list of ints, sees whether any of <s> are in <group>."""
    ints = fmt.to_list_of_ints(s)
    for val1 in ints:
        for val2 in group:
            if val1 == val2:
                return False
    return True


def internal_distance_lt(s: str, dis: int) -> bool:
    """Given that s is a list of ints, returns whether all ints are within <dis>
    of eachother."""
    ints = fmt.to_list_of_ints(s)
    for i, val1 in enumerate(ints[:-1]):
        for val2 in ints[i + 1:]:
            if abs(val1 - val2) > dis:
                return False
    return True


def is_list_of_ints(s: str) -> bool:
    """Returns whether <s> can be read as a list of integers."""
    try:
        vals = fmt.to_list_of_ints(s)
    except ValueError:
        return False
    if len(vals) == 0:
        return False
    return True


def is_valid_range(s: str, stt1: int = 0, stp1: int = 0, stt2: int = 0,
                   stp2: int = 0) -> bool:
    """Returns whether <s> can be interpreted as a range, and whether the start
    and stop indices lie in the specified ranges."""
    try:
        start, stop = fmt.to_range(s)

    except ValueError:
        log.info('User entered invalid string: "' + s + '" when asked for a '
                                                        'range.')
        return False

    if start > stop:
        log.info('User entered invalid string: "' + s + '" when asked for a '
                                                        'range.')
        return False

    if stt1 or stp1:
        if not stt1 < start < stp1:
            log.info(''.join([
                'User entered string: "', s, '" when asked for a range.',
                '\nReturning invalid because range specification were not '
                'met :', str(stt1), '<', str(start), '<', str(stp1)]))
            return False

    if stt2 or stp2:
        if not stt2 < stop < stp2:
            log.info(''.join([
                'User entered string: "', s, '" when asked for a range.',
                '\nReturning invalid because range specification were not '
                'met :', str(stt2), '<', str(stop), '<', str(stp2)]))
            return False

    return True


def is_valid_DNA(s: str) -> bool:
    """Returns whether the input sequence is a valid DNA sequence."""
    for c in s:
        if not c.upper() in 'ATCG':
            return False
    return True


def is_valid_path(s: str) -> bool:
    """Returns whether <s> is a valid filepath."""
    if s == 'DIR':
        return True
    if Path(s).exists():
        return True
    else:
        log.info(''.join(['"', s, '"', ' does not exist.']))


def valid_msa_type(s: str) -> bool:
    if s in align.KNOWN_MSA_TYPES:
        return True
    else:
        return False


def is_valid_time(s: str) -> bool:
    """Returns whether the given string can be interpreted as a valid time
    string."""
    try:
        fmt.to_time_spec(s)
    except ValueError:
        return False
    return True


VMT_S = 'Unknown MSA type. Known types include: ' + \
        ', '.join(align.KNOWN_MSA_TYPES)
VALID_MSA_TYPE = Validation(valid_msa_type, VMT_S)

VALID_INT = Validation(valid_int, 'Not a properly formatted integer.')

VALID_FLOAT = Validation(valid_float, 'Not a properly formatted float.')

VALID_RANGE = Validation(is_valid_range, "Enter a properly formatted range. "
                         + fmt.RANGE)

VALID_PATH = Validation(is_valid_path, 'This file does not exist')

VALID_DNA = Validation(is_valid_DNA, 'Input is not a valid DNA sequence.')

VALID_TIME = Validation(is_valid_time, 'Input must be formatted as a time.')
