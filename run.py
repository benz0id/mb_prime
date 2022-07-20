from collections.abc import Sequence

from config_handling.get_parameters_script import get_config_file

config_file_name = 'configs.' + get_config_file()
print(config_file_name)
config = __import__(config_file_name, fromlist=[None])


print(config.num_threads)
