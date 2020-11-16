# Empirical analysis of tunicate tree

This directory contains the data and code for the empirical analysis of the tunicate dataset of Tsagkogeorga et al. (2009) using the epistatic doublet model of Nasrallah and Huelsenbeck (2013).
These results were used to inform the choice of parameters in the simulation study.

## [`src/`](src)
Scripts for processing the alignment data and running RevBayes analyses.

Those looking to reproduce the empirical analysis in the manuscript will be most interested in `analyze_tunicates_fixed_tree.Rev`, the RevBayes analysis script to infer the posterior distribution on epistatic doublet model parameters on the RAxML tree.
Also of interest is the file `process_alignment_and_downsample.R`, which generates the 3 alignment files required to run RAxML and the fixed-tree RevBayes analysis.

Those looking to use these scripts to analyze their own datasets under the epistatic doublet model may also be interested in the file `analyze_tunicates.Rev`, which is suitable for joint inference of epistatic doublet model parameters and the phylogeny (assuming IID Exponential(10)) branch lengths.
Also of interest is the file `epistatic_doublet_model.Rev`, which contains the model for paired sites as Rev code.

## [`data/`](data)
Original alignments from Tsagkogeorga et al. (2009) and subdivisions thereof.
The file unpartitioned_50.fasta contains a subset of 50 taxa used to infer the RAxML tree.
The alignment files pair_50.nex and loop_50.nex are subdivisions of the same 50 taxa that include only paired (stem) sites or unpaired (loop) sites of the alignment, and were used to infer parameters of the epistatic doublet model.
These three alignments originate from `Tsagkogeorga-BMCEvolBiol2010-Tunicata_18S_110taxa.nex`.

## [`output`](output)
RevBayes analyses will print log files here.

## [`output/RAxML`](output/RAxML)
Maximum likelihood phylogeny from RAxML and other RAxML output.
