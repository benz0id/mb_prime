from hetero_spacer_generator.spacer_generator.multi_align_gen import *
import pytest


@pytest.mark.parametrize('bases,expected_score',
                         [
                             # ([, , , , ], ),
                             # ([, , , , ], ),
                             # ([, , , , ], ),
                             ([0, 0, 0, 0, 4], 100),
                             ([5, 5, 2, 10, 0], 256130),
                             ([0, 0, 0, 1, 1], 50 ** 2),
                             ([0, 0, 1, 1, 1], (100 / 3) ** 3),
                             ([6, 1, 2, 3, 4], 25 ** 4),
                             ([0, 1, 1, 1, 1], 25 ** 4),
                             ([4, 1, 1, 1, 1], 25 ** 4)
                         ])
def test_correct_column_score(bases, expected_score) -> None:
    """Tests that the expected score is computed from the given bases."""
    assert compute_column_score(bases, sum(bases)) == int(expected_score)


@pytest.fixture
def num_hetero() -> int:
    return 12


@pytest.mark.parametrize('align,spacers,expected_min_score', [
    ([
         'ATCTCGAATCGAT',
         'TCGAATCGAT',
         'ATCGATCGAATCGAT',
         'ATCGATCGAATCGAT',
         'TCGATCGAATCGAT',
         'CGATCGAATCGATCGA',
         'TCGATCGAATCGATCGA',
     ], [0, 1, 2, 3], 25 ** 4),
    ([
         'ATCGATCGA',
         'ATCGATCGA',
         'ATCGATCGA',
         'ATCGATCGA'
     ], [0, 1, 0, 1], 50 ** 2),
    ([
         'ATCGATCGA',
         'ATCGATCGA',
         'ATCGATCGA',
         'ATCGATCGA'
     ], [0, 1, 2, 3], 25 ** 4),
    ([
         'ATCGATCGA',
         'TCGATCGAA',
         'CGATCGATA',
         'GATCGATCA'
     ], [0, 0, 0, 0], 100),
    ([
        'ATCTGATCG',
        'TCGTATCGA',
        'CGATTCGAT',
        'GATTCGATC',
        'ATTCGATCA'
    ], [0, 0, 0, 0, 1], 100)
])
def test_correct_min(align, spacers, expected_min_score) -> None:
    """Tests that individual alignments have the expected scores."""
    lens = [len(s) for s in align]
    num_hetero = min(lens)
    aln = NumpyBindingAlign(align, spacers, num_hetero)
    assert aln.find_min_div() == expected_min_score

