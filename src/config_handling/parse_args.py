import argparse
from typing import Tuple


def get_cla_namespace() -> argparse.Namespace:
    """
    Returns a namespace storing all valid arguments given by the user.
    """

    parser = argparse.ArgumentParser()

    parser.add_argument('--config', type=str,
                        help='the config file to be used.')

    parser.add_argument('--how_random', type=int,
                        help='degree of randomness when selecting binding'
                             'sequences.')

    parser.add_argument('--num_repetitions', type=int,
                        help='number of repetitions to perform.')

    parser.add_argument('--runtime_estimate', type=int, nargs=3,
                        metavar=('HOURS', 'MINUTES', 'SECONDS'),
                        help='estimated runtime of the program.')

    parser.add_argument('--verbose', action='store_true',
                        help='whether to enable advanced output.')

    parser.add_argument('--no_warn', action='store_true',
                        help='whether to mute input-requiring warnings.')

    parser.add_argument('--silent', action='store_true',
                        help='whether to execute the program entirely.')

    parser.add_argument('--hetero_region_length', type=int,
                        help='length of heterogeneity region in primers.')

    parser.add_argument('--num_threads', type=int,
                        help='number of threads to spawn')

    parser.add_argument('--max_spacer_length', type=int,
                        help='maximum length of any spacer. Setting too short'
                             'will result in failure.')

    parser.add_argument('--out_filepath', type=str, default='',
                        help='filepath for result output')

    parser.add_argument('--max_successes', type=int, default=-1,
                        help='program will terminate after this many successful'
                             ' runs.')

    return parser.parse_args()


def override_module_variables(mod, namespace: argparse.Namespace) -> None:
    """Overrides all variables in <mod> with variables stored in <namespace>.
    <mod> must implement overwrite_var"""
    for name, value in namespace._get_kwargs():
        if value is not None:
            mod.overwrite_var(name, value)


def get_modified_config(config_module) -> None:
    """Overrides arguments in the provided config with arguments supplied by the
    user through the command line."""
    ns = get_cla_namespace()
    override_module_variables(config_module, ns)