# Overview 

This program is intended to aid in the construction of primers for metabarcoding. Given some alignment and target sites on that alignment, the 
program will select the most conserved<sup>*1</sup>, valid<sup>*2</sup> primers that amplify those sites. The program can additionally design heterogeneity spacers that
introduce diversity into the first few bases in a read, as required by some sequencing technologies. When possible, the program will attempt to reduce the propensity of primers to form dimers.


1: The conservation of nucleotides closer to the 3' end are weighted exponentially higher, due to their role in
replication initiation.

2: Primer are considered valid iff they have a GC clamp present (in the two most 3' residues), are within the melting temperature range specified by the user, and do 
not overlap with other targets on the same alignment (to avoid alternate structure formation during amplification.

# Installation

Written for python 3.10

To install the required packages, run the following.
> pip install pipreqs
> 
> pip install -r requirements.txt


# Usage
In order to start the program, execute:

>python3 run.py

This will allow you to create a config file and execute it.

# Primer Structure & Definitions

5' - Adapter - Index - SequencingPrimer - HeterogeneitySpacer - Binding Sequence - 3'

Note that the 3 most 5' components may vary depending on sequencing technology used.

Adapter - Binds the flow cell.

Index - Allows for identification of reads.

Sequencing primer - Allows for binding of the sequencing primer.

Heterogeneity spacer - Inserts diversity into the first few nucleotides in each read.

Binding sequence - Binds to regions flanking the target region in order to enable amplification of the target.

## Glossary

Target region - The portion of the target organism's genome to be amplified.

Binding Region - The regions of the target organism's genome to be bound directly by the primers.

# Parameters

The program will prompt you to specify some parameters depending on the type of run you'd like to perform.
These will be placed into a config file. Direct modification of config files is encouraged. 

## Basic Parameters

alignments_path - Path to the file containing the alignments on which the desired targets lie.

alignment_type - The filetype of the alignments in <alignments_path>.

config_type - The type of config file the user would like to create. Decides which parameters are necessary and what portions of the program to use.

adapters - The sequences to be placed 5' of the heterogeneity spacer. If using a 2-step PCR protocol, these may just contain the sequencing primer. If using a 1-step 
PCR protocol these may contain the flow cell adapter, index, and sequencing primer.

## Binding Region Selection Parameters

targets - Locations of the target sites.

target_region_len - The maximum and minimum lengths (inclusive) of any given target region.

binding_region_len - The maximum and minimum lengths (inclusive) of any given binding region.

ideal_binding_size - The ideal length of any binding region. This will serve as a start point from which adjustments will be made to reach melting temp requirements.

max_binding_target_len - The maximum combined length of the target and binding regions (forward and reverse) in any amplicon.

primer_primer_distance - The minimum distance any forward and reverse binding region. Ensures that overlap between binding regions doesn't occur.

primer_target_distance - The minimum distance between a binding region and its target.

max_mt_deviance - The maximum difference in melting temperatures between any two binding regions.

target_melting_temp - The target melting temp about, about which 

## Heterogeneity Spacer Generation Parameters

binding_sequences - The binding sequences to incorporate into the final primers, and for which heterogeneity spacers will be designed.

hetero_region_len - The number of first-read bases that must have perfect heterogeneity.

max_spacer_length - The maximum length of any heterogeneity spacer. Making this too short will often result in failure.

## Runtime Parameters

runtime_estimate - The estimated runtime of the program. The program will automatically reduce sample sizes in order to meet stringent time requirements.

num_threads - The number of child processes to spawn during any parallelizable stage of the program.

## Advanced Options

how_random - By default, the program will always choose the best option according to its own criteria. If the single best option is not entirely satisfactory, 
increasing this parameter will cause the program to randomly select one of the <how_random> best available solutions at several stages in the program.

num_repetitions - The number of times to run the program. Useful when paired with higher values of how_random generating multiple sets.

verbose - Whether to print the log file directly to console during execution.