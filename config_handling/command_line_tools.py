# funcitons for extracting and storing config files.
from typing import List
import logging
import sys
from time import sleep
import config_handling.input_validator as iv

from config_handling.input_validator import Validation

log = logging.getLogger('root')

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


def print_title(txt: str) -> None:
    """Prints <txt> formatted as a new section."""
    txt = '    ======== ' + txt + ' ========'
    print(txt)


def prompt(msg: str, validations: List[Validation]) -> str:
    """Prompts the user for input by printing <msg>, will print <err_msg> and
    retry if input does not pass <validation>."""

    valid_input_recieved = False
    inp = ''
    prompt_str = ''.join(['[?] ', msg, '\n> '])
    while not valid_input_recieved:
        valid_input_recieved = True
        inp = input(prompt_str).strip()

        # Check if user input passes
        for valid in validations:
            # Whether the valid_input_received should be changed to false.
            set_false = False

            # Check validity of input.
            is_valid = valid.is_valid(inp)
            if not is_valid:
                eprint('\n' + valid.err_msg + '\n')
                set_false = True
                sleep(0.2)

            # Prompt user to ignore.
            do_ignore_check = not is_valid and valid.is_warning
            if do_ignore_check and yes_no_prompt('Ignore?'):
                set_false = False

            if set_false:
                valid_input_recieved = False
                break

    return inp


def menu(options: List[str], title: str, instructions: str) -> int:
    """Returns the index of the item selected by the user."""
    print_title(title)
    for i, option in enumerate(options):
        print(''.join(['[', str(i + 1), '] ', option]))

    mthd = lambda s: iv.in_range_incl(s, 1, len(options))
    valid = iv.Validation(mthd, 'Invalid option number.')

    inp = prompt(instructions, [valid, iv.VALID_INT])

    return int(inp) - 1


def yes_no_prompt(msg: str) -> bool:
    """Asks the user whether they'd like to ign"""
    valid = Validation(lambda s: s.strip() in 'YN', '[Y/N]')
    s = prompt(msg + ' [Y/N]', [valid])
    if s.strip() == 'Y':
        return True
    return False
