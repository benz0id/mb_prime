import primer3
import pytest

from config_handling.formatting import AdapterPair, TargetRegionInfo
from seq_alignment_analyser.scoring import *
from seq_alignment_analyser.sequence_management import *
from test_files.fixtures_and_helpers import get_msa


msa = get_msa('alex_test_align')

@pytest.fixture
def simple_ppm() -> PrimerPartsManager:

    targ1 = TargetRegionInfo('targ1', 'alex_test_align', [50])
    targ2 = TargetRegionInfo('targ2', 'alex_test_align', [250])
    targ3 = TargetRegionInfo('targ3', 'alex_test_align', [350])
    targ4 = TargetRegionInfo('targ4', 'alex_test_align', [450])

    tt5s = [
        AdapterPair('ATATATATAATTCC', 'ATATATATAATTCC'),
        AdapterPair('ATATATATAATTCC', 'ATATATATAATTCC'),
        AdapterPair('ATATATATAATTCC', 'ATATATATAATTCC'),
        AdapterPair('ATATATATAATTCC', 'ATATATATAATTCC')
    ]

    targs = [targ1, targ2, targ3, targ4]
    msa_to_targs = {msa: targs}

    return PrimerPartsManager(tt5s, msa_to_targs)


@pytest.fixture
def simple_bp_start() -> BindingPair:
    return BindingPair(msa, 'targ1', 0, 39, 20, 20)


@pytest.fixture
def simple_bp_end() -> BindingPair:
    return BindingPair(msa, 'targ1', len(msa) - 40, len(msa) - 1, 20, 20)


@pytest.fixture
def simple_sbp(simple_ppm, simple_bp_start) -> ScoreBindingPair:
    sbp = ScoreBindingPair(simple_ppm, 65, 5)
    sbp.set_mts(simple_bp_start)
    sbp.add_bp(simple_bp_start)
    simple_ppm.add_bp(simple_bp_start)
    return sbp


@pytest.fixture
def simple_bps() -> Tuple[BindingPair, ...]:
    bp1 = BindingPair(msa, 'targ1', 102, 202, 20, 20)
    bp2 = BindingPair(msa, 'targ2', 205, 305, 20, 20)
    bp3 = BindingPair(msa, 'targ3', 310, 410, 20, 20)
    bp4 = BindingPair(msa, 'targ4', 420, 520, 20, 20)
    return bp1, bp2, bp3, bp4


def test_bp_seq(simple_bp_start, simple_bp_end) -> None:
    consensus = msa.get_consensus()
    assert simple_bp_start.f_seq + rev_comp(simple_bp_start.r_seq) == \
           consensus[:40]

    consensus = msa.get_consensus()
    assert simple_bp_end.f_seq + rev_comp(simple_bp_end.r_seq) == \
           consensus[-40:]


def test_basic_mt_functionality(simple_sbp: ScoreBindingPair, simple_bps) -> None:
    bp1, bp2, bp3, bp4, *a = simple_bps
    for bp in simple_bps:
        simple_sbp.set_mts(bp)
        seqs = simple_sbp._primer_parts.get_all_binding_seqs()
        mts = [primer3.calcTm(seq) for seq in seqs]
        mts.extend([bp.f_mt, bp1.r_mt])
        assert len(mts) == 4
        dev = simple_sbp.get_mt_deviance(bp)
        for mt in mts:
            assert abs(mt - bp.f_mt) <= abs(dev)
            assert abs(mt - bp.r_mt) <= abs(dev)


def test_eval_binding_dimer_functionality(simple_ppm, simple_bps,
                                          simple_sbp) -> None:
    bp1, *c = simple_bps
    dimer_scores = simple_sbp.get_bp_dgs(bp1)
    assert len(dimer_scores) == 7
    bp1.f_seq = 'C' * 20
    bp1.r_seq = 'G' * 20
    dimer_scores = simple_sbp.get_bp_dgs(bp1)
    last = dimer_scores[-1]
    for score in dimer_scores:
        assert score >= last


def test_eval_adapter_dimer_functionality(simple_ppm, simple_bps,
                                          simple_sbp) -> None:
    bp1, *c = simple_bps
    dimer_scores = simple_sbp.get_5p_seqs_dgs(bp1)
    assert len(dimer_scores) == 16
    simple_ppm._fivep_seqs[-1] = AdapterPair(simple_ppm._fivep_seqs[-1][0],
                                             'C' * 20)
    bp1.r_seq = 'G' * 20
    dimer_scores = simple_sbp.get_5p_seqs_dgs(bp1)
    last = dimer_scores[-1]
    for score in dimer_scores:
        assert score >= last

def test_average_conservation(simple_ppm, simple_sbp,
                              simple_bp_start, simple_bp_end) -> None:
    simple_sbp.get_avg_conservation(simple_bp_start)
    simple_sbp.get_avg_conservation(simple_bp_end)

