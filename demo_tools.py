from primer_tools import MBPrimerBuilder
from typing import List, Union, Tuple
from defaults import ABSOLUTE_MAX_NUM_HETERO, ABSOLUTE_MAX_SPACER_LENGTH
from presenters import Presenter
from time import sleep

STR = 'str'
RANGE = 'ran'
DNA = 'dna'


def valid_input(input: Union[int, str], mode: str, start: Union[int, str] = 0,
                end: Union[int, str] = 0,
                allowed: Union[int, str] = '') -> bool:
    """Returns whether the input is an allowed input according to some
    parameters.
    Modes:
    STR:
        <input> is valid iff <input> in <allowed>
    RANGE:
        <input> is valid iff <start> =< <input> <= <end>.
        all of <input>, <start>, and <end> must be integers or string integers.
        """
    if mode == RANGE:
        start = int(start)
        end = int(end)
        return start <= int(input) <= end
    if mode == STR:
        return str(input) in str(allowed)
    if mode == DNA:
        for char in input:
            if char.capitalize() not in 'ATCG':
                return False
        return True


def while_not_valid(msg: str, e_msg: str, mode: Union[int, str],
                    start: Union[int, str] = '0',
                    end: Union[int, str] = '0',
                    allowed: Union[int, str] = '') -> str:
    """Continually prompts the user to enter input with the given <e_msg> while
    valid_input run with the given parameters is not true. Will try once with
    <msg> before printing <e_msg>."""
    inpt = input(msg)
    while not valid_input(inpt, mode, start, end, allowed):
        inpt = input(e_msg)
    return inpt


def get_incomplete_primer() -> MBPrimerBuilder:
    """Prompts the user through the console for the components required to build
    a metabarcoding primer."""
    primer = MBPrimerBuilder()
    inpt = while_not_valid(
        "Enter the adapter region of your primer (5' - 3'): ",
        "Bad input, please ensure your entry is a valid sequence:", DNA)
    primer.set_adapter_seq(inpt.upper())

    inpt = while_not_valid(
        "Enter the indexing region of your primer (5' - 3'): ",
        "Bad input, please ensure your entry is a valid sequence:", DNA)
    primer.set_index_seq(inpt.upper())

    inpt = while_not_valid(
        "Enter the binding region of your primer (5' - 3'): ",
        "Bad input, please ensure your entry is a valid sequence:", DNA)
    primer.set_binding_seq(inpt.upper())

    return primer


def get_params() -> Tuple[int, int]:
    """Prompts the user through the console for the parameters required to build
    a metabarcoding primer.
    Returns: Tuple[max spacer length, heterogeneity region size]"""
    spacer_length = while_not_valid(
        "Enter the maximum size of any heterogeneity spacer to be generated: ",
        "Bad input, please ensure your entry is a number greater than 0:",
        RANGE, start=0, end=ABSOLUTE_MAX_SPACER_LENGTH)

    num_hetero = while_not_valid(
        "Enter the size of the heterogeneity region you'd like in your primer: ",
        "Bad input, please ensure your entry is a number greater than 0:",
        RANGE, start=0, end=ABSOLUTE_MAX_NUM_HETERO)

    return int(spacer_length), int(num_hetero)


def print_dots(delay: int) -> None:
    """Prints a dot every <delay> seconds."""
    while True:
        sleep(delay)
        print('.', end='')

