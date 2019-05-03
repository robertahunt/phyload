# phyload simulation study
Can epistasis break phylogenetic models? If so, where?

## Dependencies
- [Conda](https://conda.io/):
After installing Conda, we recommend setting up an environment named _phyload_ with the required version of Python as follows:
```bash
conda create -n phyload python=3.7
conda activate phyload
```
- [Biopython](https://biopython.org)
```bash
conda install -c conda-forge biopython
```
- [SCons](https://scons.org)
```bash
conda install scons
```
- [RevBayes](https://revbayes.github.io/software)
Follow the install instructions.

## SCons pipeline: [SConstruct](SConstruct)
Perform the simulation with default parameters by running
```bash
$ scons
```
For help on command line arguments, run
```bash
$ scons -h
```
and note `Local Options` at the bottom.

## [`simulation_scripts/`](simulation_scripts)
Scripts needed to simulate a single cell in the simulation study.

### [`mk_alns.R`](simulation_scripts/mk_alns.R)

To perform a single simulation, use
```bash
$ Rscript simulation_scripts/mk_alns.R <n_iid> <n_epi> <d> <seed> <outbase> <revpath>
```
Separate iid and epistatic alignments will be written in nexus format to `<outbase>/iid_aln.nex` and `<outbase>/epi_aln.nex`.

The arguments are
- `n_iid`: number of iid sites
- `n_epi`: number of epistatically paired sites (must be even)
- `d`: *d* parameter for Nasrallah-Huelsenbeck model
- `seed`: random seed passed to RevBayes
- `outbase`: path to write output alignments
- `revpath`: path to installed version of RevBayes (can be `rb` if it is on `$PATH`)

### Utilities

- [`merge_alns.R`](simulation_scripts/merge_alns.R)
Merge an iid nexus alignment `<iid_aln>` and and epistatic nexus alignment `<epi_aln>` (as produced by `mk_alns.R`) and write to a new nexus alignment `<merged_aln>` with
```bash
$ Rscript simulation_scripts/merge_alns.R <iid_aln> <epi_aln> <merged_aln>
```
- [`epistatic_doublet_model_stub.Rev`](simulation_scripts/epistatic_doublet_model_stub.Rev) Core Rev script for the epistatic model.
- [`rev_model_template.Rev`](simulation_scripts/rev_model_template.Rev) A template for a Rev script to simulate 100 alignments, missing key parameters.
- [`simulation_tree.tre`](simulation_scripts/simulation_tree.tre) The treefile used for simulating alignments.





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
