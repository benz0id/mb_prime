import collections
import logging
from typing import List, Tuple

log = logging.getLogger('root')

InclRange = collections.namedtuple('InclRange', ['start', 'stop'])

AdapterPair = collections.namedtuple('AdapterPair', ['forward', 'reverse'])

TimeSpec = collections.namedtuple("TimeSpec", ['hours', 'seconds', 'minutes'])

Sites = Tuple[int, ...]

TargetRegionInfo = collections.namedtuple('TargetRegionInfo',
                                          ['name', 'aln_filename', 'sites'])

PrimerParams = collections.namedtuple(
    'PrimerParams',
    [
        'primer_primer_distance', 'primer_target_distance',
        'target_region_len', 'binding_region_len', 'ideal_binding_size',
        'max_binding_target_len'
    ])

INT_LIST = 'Format: "val1, val2, val3..."'
RANGE = 'Format: "start, stop" (inclusive)'
TIME = 'Format "int: int: int" (hrs:mins:secs)'


def incl_to_range(incl_range: InclRange) -> range:
    """Converts the given <incl_range> to an actual range."""
    return range(incl_range.start, incl_range.stop + 1)


def int_strip(s: str) -> int:
    """Removes all non-numeric characters from s and converts it to an int.
    Raises ValueError upon failure. Integer must be continuous series of
    numbers."""
    int_s = ''
    int_start = False
    int_fin = False

    for c in s:
        if c.isnumeric() or c in '-':
            if int_start and int_fin:
                log.info('User entered string: "' + s +
                         '" when asked for an int.')
                raise ValueError('Integer is not continuous.')
            if not int_start:
                int_start = True

            int_s += c
        else:
            if int_start and not int_fin:
                int_fin = True
    if not int_s:
        log.info('User entered string: "' + s + '" when asked for an int.')
        raise ValueError('String did not contain an int.')

    if '-' in s and s[0] != '-':
        log.info('User entered string: "' + s + '" when asked for an int.')
        raise ValueError('String did not contain an int.')

    return int(int_s)


def to_list_of_ints(s: str) -> List[int]:
    """Converts the given string to a list of integers if possible. Raises
    ValueError otherwise."""
    vals = []
    int_strs = s.split(',')
    for int_str in int_strs:
        vals.append(int_strip(int_str))
    return vals


def to_range(s: str) -> Tuple[int, int]:
    """Convers the given string to a range. Raises ValueError on failure."""
    err = ValueError('Inappropriate range format given.')

    vals = s.split(',')
    if len(vals) != 2:
        raise err

    try:
        vals = int_strip(vals[0]), int_strip(vals[1])
    except ValueError:
        raise err

    return vals


def to_time_spec(s: str) -> TimeSpec:
    """Converts the given string to a timespec."""
    strs = s.split(':')
    if len(strs) != 3:
        raise ValueError('Incorrect formatting.')
    vals = []
    for val in strs:
        vals.append(int_strip(val))

    return TimeSpec(vals[0], vals[1], vals[2])








