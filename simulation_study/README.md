# phyload simulation study
Can epistasis break phylogenetic models? If so, where?

## [`simulation_scripts`](simulation_scripts)
Scripts needed to simulate a single cell in the simulation study, and a master script to use them, `simulate_grid_cell.R`.

To simulate a cell, from top level of phyload repository, use `Rscript simulate_grid_cell.R nsites prop_epi d "relative/output/path" "path/to/rb" seed`.

The arguments are
- nsites: number of total sites
- prop_epi: proportion of nsites that are part of paired epistatic interactions
- d: d parameter for Nasrallah-Huelsenbeck model
- "relative/output/path": path to output alignments (relative to location of phyload repository)
- "path/to/rb": absolute file path to installed version of RevBayes (can be "rb" if it is installed)
- seed: seed passed to RevBayes

The other files in this folder are required for the above to work, and any change to file names will cause the setup to break.
- [`simulation_scripts/epistatic_doublet_model_stub.Rev`](simulation_scripts/epistatic_doublet_model_stub.Rev) Core Rev script for the epistatic model.
- [`simulation_scripts/merge_alignments.R`](simulation_scripts/merge_alignments.R) The actual simulation produces two alignments, one with (1 - prop_epi) sites from a plain GTR model, one with 0.5 * prop_epi sites, each of which is a pair of interacting sites. This script unfolds the pairs into sites and makes a single alignment out of everything
- [`simulation_scripts/rev_model_template.Rev`](simulation_scripts/rev_model_template.Rev) A template for a Rev script to simulate 100 alignments, missing key parameters.
- [`simulation_scripts/simulation_tree.tre`](simulation_scripts/simulation_tree.tre) The treefile used for simulating alignments.

## [`preliminary`](preliminary)
This is for our rough, first-pass simulation attack

Here we simulate alignments from Nasrallah and Huelsenbeck's 2013 model of length _n_.
We simulate a 4x5 grid, with a proportion of sites (0.25,0.5,0.75,1.0) being drawn from the epistasis model with d parameter (0,1,2,5,10).
When d is 0, there is no real epistatic interactions and it is just a model with some odd pairwise stationary frequencies.
As our references, we simulate alignments of length (1,0.75,0.5,0.25) x _n_ as a reference.
The question of epistasis breaking the model then can be addressed by comparing the accuracy of trees estimated from epistatic alignments to the accuracy of the trees estimated from these purely site-IID models.
We take _n_ = 1496 because it's about the size of flu HA alignments and is divisible by 8 (if a proportion _p_ of the alignment is drawn from the epistatic model, there must be 0.5_pn_ sites, each paired with the other 0.5_pn_ sites.)

### [`preliminary/1_decreasing_site_counts`](preliminary/1_decreasing_site_counts)

Contains the reference alignments of decreasing numbers of sites. Subdirectories named frac-<_x_\> where _x_ is the fraction of the full alignment size _n_ (all alignments in this directory are this size).

### [`preliminary/2_epistasis`](preliminary/2_epistasis)

Contains the alignments with epistatic portions.
Subdirectories named prop-<_x_\>-d-<_y_\>, where _x_ is the proportion of sites simulated form the epistatic interaction model with parameter _d_ = _y_ (all alignments in this directory are simulated according to this model).

### [`preliminary/0_setup`](preliminary/0_setup)
Directory with files and scripts for simulating and analyzing results.
Currently all scripts live in the `src` subdirectory.
Within this,
  - `analysis_scripts` will contain Rev scripts to analyze the simulated data
  - `generate_analysis_scripts` will contain scripts to generate the analysis scripts
  - `simulation` contains the scripts required to simulate the alignments
  - `slurm_scripts` will contain scripts to automate running the analyses
  - `summarization` will contain scripts to process the analysis results, diagnose MCMC performance, etc.

When simulating alignments using `simulate_all_alignments.R`, it is recommended that simulate_all_alignments.R is used with Rscript, it can hang in Rstudio.
