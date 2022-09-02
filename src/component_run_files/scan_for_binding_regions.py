from src.execution_managers.demo_tools import while_not_valid, RANGE
from src.execution_managers.parameter_manager import get_pm
from src.hetero_spacer_generator.defaults import NUM_TO_KEEP
from seq_alignment_analyser.align import MSA
from seq_alignment_analyser import BestPrimers, get_region, \
    get_seqs, vis_score


def main():
    pm = get_pm()

    # Ensure that necessary parameters have been defined.
    if not pm.can_sfbr():
        raise ValueError('Ensure all required parameters in the config file are ' +
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
    GRAPH_OUT_FILEPATH = pm.get('GRAPH_OUT_FILEPATH')
    STORE_IMAGE = pm.get('STORE_IMAGE')
    DISPLAY_GRAPH = pm.get('DISPLAY_GRAPH')
    WINDOW_SIZE = pm.get('WINDOW_SIZE')

    # Some error checking on inputs.
    cond1 = FORWARD_BINDING_START <= FORWARD_BINDING_END
    cond2 = REVERSE_BINDING_START <= REVERSE_BINDING_END
    cond3 = AMPLICON_LENGTH_MIN <= AMPLICON_LENGTH_MAX

    if not (cond1 and cond2 and cond3):
        raise ValueError('Ensure all specified ranges are valid.')

    print('Close Graph to Begin Analysis.')

    # Extract specified ranges and construct Alignment.
    alignment = MSA(MSA_FILEPATH, WINDOW_SIZE)
    f_allowed_5p = list(range(FORWARD_BINDING_START, FORWARD_BINDING_END + 1))
    r_allowed_5p = list(range(REVERSE_BINDING_START, REVERSE_BINDING_END + 1))
    a_allowed_len = list(range(AMPLICON_LENGTH_MIN, AMPLICON_LENGTH_MAX))

    # Display alignment graph.
    alignment.gen_plot()
    alignment.add_primers_to_graph(FORWARD_BINDING_START, FORWARD_BINDING_END,
                                   REVERSE_BINDING_START, REVERSE_BINDING_END,
                                   primer_region_name="Potential Primer Region",
                                   amplicon_region_name='Guaranteed Amplicon')
    if DISPLAY_GRAPH:
        alignment.show_plot()

    if STORE_IMAGE:
        alignment.save_plot(GRAPH_OUT_FILEPATH)

    # Perform analysis.
    bp = BestPrimers(alignment, f_allowed_5p, r_allowed_5p, FORWARD_BINDING_LENGTHS,
                     REVERSE_BINDING_LENGTHS, a_allowed_len, FORWARD_ADAPTERS,
                     REVERSE_ADAPTERS, NUM_TO_KEEP)

    best = bp.get_n_best(5, V=True)

    print("Analysis complete...")

    for i, bps in enumerate(best):
        print('BINDING PAIR', i + 1)
        print(vis_score(bps))

    # Allow user to anaylse selected options.
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
        f_start, f_stop = get_region(fp, False)
        print(alignment.scan_region(f_start, f_stop + 1))

        print('\nReverse Binding Region Analysis:')
        r_start, r_stop = get_region(rp, False, 'r')
        print(alignment.scan_region(r_start, r_stop + 1))

        # Display alignment graph.
        alignment.gen_plot()
        alignment.add_primers_to_graph(f_start, f_stop, r_start, r_stop)
        if DISPLAY_GRAPH:
            alignment.show_plot()

        if STORE_IMAGE:
            alignment.save_plot(GRAPH_OUT_FILEPATH)

        if __name__ != '__main__':
            break

if __name__ == "__main__":
    main()












