# Stores critera used by the rest of the program to evaluate potential primers.
from typing import Dict, List, Tuple
import abc
from Bio.Seq import Seq
from hetero_spacer_generator.defaults import CONSEC_WEIGHT, TOTAL_WEIGHT
from hetero_spacer_generator.primer_tools import HeteroSeqTool, \
    MBPrimer, MBPrimerBuilder, PairWiseCriterionSingle, rev_seq
from hetero_spacer_generator.sequence_tools import SeqAnalyzer, get_p3_adapter


class EvalMBPrimer(HeteroSeqTool):
    """Generates criteria to score primers on their propensity to form secondary
    structures."""

    def __init__(self, max_spacer_length: int, num_hetero: int):
        HeteroSeqTool.__init__(self, max_spacer_length, num_hetero)

    @abc.abstractmethod
    def get_homo_criteria(self) -> \
            Tuple[List[PairWiseCriterionSingle],
                  Dict[PairWiseCriterionSingle, int]]:
        """Returns a tuple containing a list and a dictionary. The list contains
        callables that score *individual* metabarcoding primers on their ability to
        form homodimers, focusing on the effects of the heterogenity region.
        The dictionary maps each of these criteria to their weights. primer_eval is
        passed so that users may specify additional parameters within it.
        """
        pass

    @abc.abstractmethod
    def get_hetero_criteria(self) -> Tuple[List[PairWiseCriterionSingle],
              Dict[PairWiseCriterionSingle, int]]:
        """Returns a tuple containing a list and a dictionary. The list contains
        callables that score *pairs* of metabarcoding primers on their ability to
        form homodimers, focusing on the effects of the heterogenity region.
        The dictionary maps each of these criteria to their weights. primer_eval is
        passed so that users may specify additional parameters within it.
        """
        pass


class EvalMBPrimer3(EvalMBPrimer):
    """Generates primer evaluation criteria using thermal estimates produced by
     Primer3"""

    p3a = get_p3_adapter()

    def __init__(self, max_spacer_length: int, num_hetero: int,
                 degen: bool = None):
        EvalMBPrimer.__init__(self, max_spacer_length, num_hetero)


    def start_counting_comparisons(self) -> None:
        """Begins counting the number of sequence analyses performed."""
        self.p3a.start_counting_comparisons()

    def stop_counting_comparisons(self) -> int:
        """Returns the number of sequence analyses performed since
        start_counting_comparisons was last called."""
        return self.p3a.stop_counting_comparisons()

    def get_homo_criteria(self) -> \
            Tuple[List[PairWiseCriterionSingle],
                  Dict[PairWiseCriterionSingle, int]]:
        """Returns a tuple containing a list and a dictionary. The list contains
        callables that score *individual* metabarcoding primers on their ability to
        form homodimers, focusing on the effects of the heterogenity region.
        The dictionary maps each of these criteria to their weights. primer_eval is
        passed so that users may specify additional parameters within it.
        """
        criteria = [
            self.p3a.calc_hairpin_score,
            self.p3a.calc_homodimer_score
        ]
        # Weights somewhat arbitrary.
        criteria_weights = {
            self.p3a.calc_hairpin_score: 1,
            self.p3a.calc_homodimer_score: 4
        }

        return criteria, criteria_weights

    def get_hetero_criteria(self) -> Tuple[List[PairWiseCriterionSingle],
                                           Dict[PairWiseCriterionSingle, int]]:
        """Returns a tuple containing a list and a dictionary. The list contains
        callables that score *pairs* of metabarcoding primers on their ability to
        form homodimers, focusing on the effects of the heterogenity region.
        The dictionary maps each of these criteria to their weights. primer_eval is
        passed so that users may specify additional parameters within it.
        """
        criteria = [
            self.p3a.calc_heterodimer_score,
        ]
        # Weights somewhat arbitrary.
        criteria_weights = {
            self.p3a.calc_heterodimer_score: 1,
        }

        return criteria, criteria_weights


class EvalMBPrimerNaive(EvalMBPrimer, SeqAnalyzer):
    """A class designed to evaluate the characteristics of MBPrimers"""

    _forward_primer: MBPrimerBuilder
    _rev_primer: MBPrimerBuilder

    def __init__(self, max_spacer_length: int, num_hetero: int,
                 degen: bool = None):
        EvalMBPrimer.__init__(self, max_spacer_length, num_hetero)
        SeqAnalyzer.__init__(self, degen)

    def get_homo_criteria(self) -> \
            Tuple[List[PairWiseCriterionSingle],
                  Dict[PairWiseCriterionSingle, int]]:
        """Returns a tuple containing a list and a dictionary. The list contains
        callables that score *individual* metabarcoding primers on their ability to
        form homodimers, focusing on the effects of the heterogenity region.
        The dictionary maps each of these criteria to their weights. primer_eval is
        passed so that users may specify additional parameters within it.
        """
        criteria = [
                # The number of consecutive complementary bases between the
                # heterogeneity spacer and any other region.
                self.eval_homo_hetero_spacer_binding_consec,
                # The total number of complementary bases.
                self.eval_homo_hetero_spacer_binding_total
            ]
        # Weights somewhat arbitrary.
        criteria_weights = {
            # The number of consecutive complementary bases is by defenition less
            # than or equal to the total complementarity. It will often be far
            # less, so account for this by giving it a higher weight.
            self.eval_homo_hetero_spacer_binding_consec: CONSEC_WEIGHT,
            self.eval_homo_hetero_spacer_binding_total: TOTAL_WEIGHT
        }

        return criteria, criteria_weights

    def get_hetero_criteria(self) -> Tuple[List[PairWiseCriterionSingle],
                                           Dict[PairWiseCriterionSingle, int]]:
        """Returns a tuple containing a list and a dictionary. The list contains
        callables that score *pairs* of metabarcoding primers on their ability to
        form homodimers, focusing on the effects of the heterogenity region.
        The dictionary maps each of these criteria to their weights. primer_eval is
        passed so that users may specify additional parameters within it.
        """
        criteria = [
                # The number of consecutive complementary bases between the
                # heterogeneity spacer and any other region.
                self.eval_hetero_hetero_spacer_binding_consec,
                # The total number of complementary bases.
                self.eval_hetero_hetero_spacer_binding_total
            ]
        # Weights somewhat arbitrary.
        criteria_weights = {
            # The number of consecutive complementary bases is by definition less
            # than or equal to the total complementarity. It will often be far
            # less, so account for this by giving it a higher weight.
            self.eval_hetero_hetero_spacer_binding_consec: CONSEC_WEIGHT,
            self.eval_hetero_hetero_spacer_binding_total: TOTAL_WEIGHT
        }

        return criteria, criteria_weights

    def eval_inherent_heterodimer_consec(self, forward_primer: MBPrimerBuilder,
                                         reverse_primer: MBPrimerBuilder) -> int:
        """Returns the greatest number of consecutively complementary bases
        between the components of the <forward_primer> and <reverse_primer>
        external to the heterogeneity spacer. All primers should be input 5'-3'.
        """
        forward_regions = [forward_primer.get_5p(), forward_primer.get_3p()]
        reverse_regions = [reverse_primer.get_5p(), reverse_primer.get_3p()]
        comp_method = self.get_consec_complementarity
        scores = []
        for for_region in forward_regions:
            for rev_region in reverse_regions:
                scores.append(
                    self.comp_seqs_any_overlap(for_region, rev_region,
                                               comp_method))
        return max(scores)

    def eval_inherent_homodimer_consec(self, primer: MBPrimerBuilder) -> int:
        """Returns the greatest number of consecutively complementary bases between
        the components of the <forward_primer> external to the heterogeneity
        spacer.All primers should be input 5'-3'"""
        return self.eval_inherent_heterodimer_consec(primer, primer)

    def eval_homo_hetero_spacer_binding_consec(self, primer: MBPrimer) -> int:
        """Evaluates the ability of <primer>'s heterogeneity spacer to form
        consecutive bonds with itself. All primers should be input 5'-3'"""
        f_skip = r_skip = min(len(primer.get_5p()), len(primer.get_3p()))
        return self.comp_seqs_any_overlap(primer, primer,
                                          self.get_consec_complementarity,
                                          f_skip, r_skip)

    def eval_hetero_hetero_spacer_binding_consec(self, for_primer: MBPrimer,
                                                 rev_primer: MBPrimer) -> int:
        """Evaluates the ability of <for_primer> and <rev_primer> heterogeneity spacer to form
       consecutive bonds with itself. All primers should be input 5'-3'"""
        f_skip = min(for_primer.get_5p_len(), rev_primer.get_5p_len())
        r_skip = min(for_primer.get_3p_len(), rev_primer.get_3p_len())
        return self.comp_seqs_any_overlap(for_primer, rev_primer,
                                          self.get_consec_complementarity,
                                          f_skip, r_skip)

    def eval_inherent_heterodimer_total(self, forward_primer: MBPrimerBuilder,
                                        reverse_primer: MBPrimerBuilder) -> int:
        """Returns the greatest number consecutively complementary bases between
        the components of the <forward_primer> and <reverse_primer> external to
        the heterogeneity spacer. All primers should be input 5'-3'"""
        forward_regions = [forward_primer.get_5p(), forward_primer.get_3p()]

        # Change reverse regions to 3' - 5' to allow for comparison.
        reverse_regions = [rev_seq(reverse_primer.get_5p()),
                           rev_seq(reverse_primer.get_3p())]
        comp_method = self.get_non_consec_complementarity
        scores = []
        for for_region in forward_regions:
            for rev_region in reverse_regions:
                scores.append(
                    self.comp_seqs_any_overlap(for_region, rev_region,
                                               comp_method))
        return max(scores)

    def eval_inherent_homodimer_total(self, primer: MBPrimerBuilder) -> int:
        """Returns the greatest number complementary bases between
        the components of the <forward_primer> external to the heterogeneity
        spacer. All primers should be input 5'-3'"""
        return self.eval_inherent_heterodimer_consec(primer, primer)

    def eval_homo_hetero_spacer_binding_total(self, primer: MBPrimer) -> int:
        """Evaluates the ability of <primer>'s heterogeneity spacer to form
        bonds with itself."""
        f_skip = r_skip = min(len(primer.get_5p()), len(primer.get_3p()))
        return self.comp_seqs_any_overlap(primer, primer,
                                          self.get_non_consec_complementarity,
                                          f_skip, r_skip)

    def eval_hetero_hetero_spacer_binding_total(self, for_primer: MBPrimer,
                                                rev_primer: MBPrimer) -> int:
        """Evaluates the ability of <primer>'s heterogeneity spacer to form
        bonds with itself."""
        f_skip = min(for_primer.get_5p_len(), rev_primer.get_5p_len())
        r_skip = min(for_primer.get_3p_len(), rev_primer.get_3p_len())
        return self.comp_seqs_any_overlap(for_primer, rev_primer,
                                          self.get_non_consec_complementarity,
                                          f_skip, r_skip)

