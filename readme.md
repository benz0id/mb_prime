
Complete Primer Structure:

5' - Adapter - Index - SequencingPrimer - HeterogeneitySpacer - Binding Sequence - 3'

Note that the 3 most 5' components may vary depending on sequencing technology used.

Adapter - Binds the flow cell.

Index - Allows for identification of reads.

Sequencing primer - Allows for binding of the sequencing primer.

Heterogeneity spacer - Inserts diversity into the first few nucleotides in each read.

Binding sequence - Binds to regions flanking the target region in order to enable amplification of the target.

Glossary:

Target region - The portion of the target organism's genome to be amplified.
Binding Region - The regions of the target organism's genome to be bound directly by the primers.


Basic Parameters:

alignments_path - Path to the file containing the alignments on which the desired targets lie.

alignment_type - The filetype of the alignments in <alignments_path>.

Binding Region Selection Parameters

target_region_len - The maximum and minimum lengths (inclusive) of any given target region.

binding_region_len - The maximum and minimum lengths (inclusive) of any given binding region.

ideal_binding_size - The ideal length of any binding region. This will serve as a start point from which adjustments will be made to reach melting temp requirements.

max_binding_target_len - The maximum combined length of the target and binding regions (forward and reverse) in any amplicon.

primer_primer_distance - The minimum distance any forward and reverse binding region. Ensures that overlap between binding regions doesn't occur.

primer_target_distance - The minimum distance between a binding region and its target.

max_mt_deviance - The maximum difference in melting temperatures between any two binding regions.

target_melting_temp - The target melting temp about, about which 

hetero_region_len - The number of first-read bases that must have perfect heterogeneity.

max_spacer_length - The maximum length of

runtime_estimate = TimeSpec(hours=0, minutes=2, seconds=0)

num_threads = 8

Usage:

