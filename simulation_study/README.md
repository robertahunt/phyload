# phyload simulation study
Can epistasis break phylogenetic models? If so, where?

## Dependencies
- [Conda](https://conda.io/):
After installing Conda, we recommend setting up an environment named _phyload_ with the required version of Python and R as follows:
```bash
conda create -n phyload python=3.7 r-essentials r-base
conda activate phyload
```
- [SCons](https://scons.org), [biopython](https://biopython.org), and [seaborn](https://seaborn.pydata.org)
```bash
conda install scons biopython seaborn
```
- [Nestly](https://nestly.readthedocs.io/en/latest/)
```bash
pip install nestly
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
