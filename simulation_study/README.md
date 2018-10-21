# phyload simulation study
Can epistasis break phylogenetic models? If so, where?

## [`preliminary`](preliminary)
This is for our rough, first-pass simulation attack

Here we simulate alignments from Nasrallah and Huelsenbeck's 2013 model of length _n_.
We simulate a 4x5 grid, with a proportion of sites (0.25,0.5,0.75,1.0) being drawn from the epistasis model with d parameter (0,1,2,5,10).
When d is 0, there is no real epistatic interactions and it is just a model with some odd pairwise stationary frequencies.
As our references, we simulate alignments of length (1,0.75,0.5,0.25) x _n_ as a reference.
The question of epistasis breaking the model then can be addressed by comparing the accuracy of trees estimated from epistatic alignments to the accuracy of the trees estimated from these purely site-IID models.
We take _n_ = 1496 because it's about the size of flu HA alignments and is divisible by 8 (if a proportion _p_ of the alignment is frawn from the epistatic model, there must be 0.5_pn_ sites, each paired with the other 0.5_pn_ sites.)

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
