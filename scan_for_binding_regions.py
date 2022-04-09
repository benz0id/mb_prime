from execution_managers.demo_tools import while_not_valid, RANGE
from execution_managers.parameter_manager import get_pm
from hetero_spacer_generator.defaults import NUM_TO_KEEP
from seq_alignment_analyser.align import MSA
from seq_alignment_analyser.best_primers import BestPrimers, get_region, \
    get_seqs, vis_score

pm = get_pm()

# Ensure that necessary parameters have been defined.
if not pm.can_sfbr():
    raise ValueError('Ensure all required parameters in the config file are '
                     'properly defined. Currently missing values include:' +
                     str(pm.missing_sfbr())[1:-1])

# Unpack paramaters
MSA_FILEPATH = pm.get('MSA_FILEPATH')
FORWARD_ADAPTERS = pm.get('FORWARD_ADAPTERS')
REVERSE_ADAPTERS = pm.get('REVERSE_ADAPTERS')
FORWARD_BINDING_START = pm.get('FORWARD_BINDING_START')
FORWARD_BINDING_END = pm.get('FORWARD_BINDING_END')
REVERSE_BINDING_START = pm.get('REVERSE_BINDING_START')
REVERSE_BINDING_END = pm.get('REVERSE_BINDING_END')
FORWARD_BINDING_LENGTHS = pm.get('FORWARD_BINDING_LENGTHS')
REVERSE_BINDING_LENGTHS = pm.get('REVERSE_BINDING_LENGTHS')
AMPLICON_LENGTH_MIN = pm.get('AMPLICON_LENGTH_MIN')
AMPLICON_LENGTH_MAX = pm.get('AMPLICON_LENGTH_MAX')
CONSERVATION_WEIGHT = pm.get('CONSERVATION_WEIGHT')
DIMER_WEIGHT = pm.get('DIMER_WEIGHT')

# Some error checking on inputs.
cond1 = FORWARD_BINDING_START <= FORWARD_BINDING_END
cond2 = REVERSE_BINDING_START <= REVERSE_BINDING_END
cond3 = AMPLICON_LENGTH_MIN <= AMPLICON_LENGTH_MAX

if not (cond1 and cond2 and cond3):
    raise ValueError('Ensure all specified ranges are valid.')

alignment = MSA(MSA_FILEPATH)
f_allowed_5p = list(range(FORWARD_BINDING_START, FORWARD_BINDING_END + 1))
r_allowed_5p = list(range(REVERSE_BINDING_START, REVERSE_BINDING_END + 1))
a_allowed_len = list(range(AMPLICON_LENGTH_MIN, AMPLICON_LENGTH_MAX))

bp = BestPrimers(alignment, f_allowed_5p, r_allowed_5p, FORWARD_BINDING_LENGTHS,
                 REVERSE_BINDING_LENGTHS, a_allowed_len, FORWARD_ADAPTERS,
                 REVERSE_ADAPTERS, NUM_TO_KEEP)

best = bp.get_n_best(5)

for i, bps in enumerate(best):
    print('BINDING PAIR', i + 1)
    print(vis_score(bps))

while True:
    sel = while_not_valid('Enter parameter number to extract binding sequences'
                          ' and for brief analysis.',
                          'Ensure that the given input is valid.',
                          mode=RANGE, start=1, end=5)
    sel_params = best[int(sel) - 1]
    print(get_seqs(sel_params, consensus=alignment.get_consensus()))

    fp = sel_params.f_params
    rp = sel_params.r_params
    print('\nForward Binding Region Analysis:')
    start, stop = get_region(fp, False)
    print(alignment.scan_region(start, stop + 1))

    print('\nReverse Binding Region Analysis:')
    start, stop = get_region(rp, False, 'r')
    print(alignment.scan_region(start, stop + 1))












