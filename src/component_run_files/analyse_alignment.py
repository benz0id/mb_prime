from src.execution_managers.parameter_manager import get_pm
from seq_alignment_analyser.align import MSA

pm = get_pm()

# Ensure that necessary parameters have been defined.
if not pm.can_aa():
    raise ValueError('Ensure all required parameters in the config file are '
                     'properly defined. Currently missing values include:' +
                     str(pm.missing_aa())[1:-1])

# Unpack paramaters
MSA_FILEPATH = pm.get('MSA_FILEPATH')
GRAPH_OUT_FILEPATH = pm.get('GRAPH_OUT_FILEPATH')
STORE_IMAGE = pm.get('STORE_IMAGE')
DISPLAY_GRAPH = pm.get('DISPLAY_GRAPH')
WINDOW_SIZE = pm.get('WINDOW_SIZE')

msa = MSA(MSA_FILEPATH, WINDOW_SIZE)

msa.gen_plot()

if DISPLAY_GRAPH:
    msa.show_plot()

if STORE_IMAGE:
    msa.save_plot(GRAPH_OUT_FILEPATH)





