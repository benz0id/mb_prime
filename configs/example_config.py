import os
from src.config_handling.formatting import *
from pathlib import Path
import inspect
from typing import Any

def overwrite_var(name: str, val: Any) -> None:
    f = inspect.currentframe()
    f.f_globals[name] = val

# These are the parameters that I'd recommend playing with.

how_random = 3

config_type = 'binding'

num_repetitions = 5

max_spacer_length = 8

hetero_region_len = 10

runtime_estimate = TimeSpec(hours=0, minutes=10, seconds=0)

verbose = False

# =============================================================================

alignments_path = Path(os.path.dirname(__file__)).parent / 'alignments'

alignment_type = 'fasta'

target_region_len = InclRange(start=0, stop=150)

binding_region_len = InclRange(start=15, stop=20)

ideal_binding_size = 20

max_binding_target_len = 126

target_melting_temp = 55.0

max_mt_deviance = 5.0

primer_primer_distance = 0

primer_target_distance = 5

num_threads = 8

targets = [
	TargetRegionInfo(name='targ_1', aln_filename='example_alignment.fas', sites=(247, 249)),
	TargetRegionInfo(name='targ_2', aln_filename='example_alignment.fas', sites=(364, 372)),
	TargetRegionInfo(name='targ_3', aln_filename='example_alignment.fas', sites=(781, 783)),
	TargetRegionInfo(name='targ_4', aln_filename='example_alignment.fas', sites=(874, 897))
	]

slc = slice(None)

adapters = [
	SeqPairParam('ACACTCTTTCCCTACACGACGCTCTTCCGATCT'[slc], 'GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT'[slc]),
	SeqPairParam('ACACTCTTTCCCTACACGACGCTCTTCCGATCT'[slc], 'GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT'[slc]),
	SeqPairParam('ACACTCTTTCCCTACACGACGCTCTTCCGATCT'[slc], 'GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT'[slc]),
	SeqPairParam('ACACTCTTTCCCTACACGACGCTCTTCCGATCT'[slc], 'GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT'[slc])
	]

out_filepath = ''