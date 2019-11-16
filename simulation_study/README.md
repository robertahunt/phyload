# phyload simulation study
Can epistasis break phylogenetic models? If so, where?

## Dependencies
- [Conda](https://conda.io/):
After installing Conda, we recommend setting up an environment named _phyload_ with the required versions of R and Python, and packages [SCons](https://scons.org), [biopython](https://biopython.org), [seaborn](https://seaborn.pydata.org), [phangorn](https://cran.r-project.org/web/packages/phangorn/index.html), and [coda](https://cran.r-project.org/web/packages/coda/coda.pdf):  
```bash
conda create -n phyload python=3.7 r-essentials r-base scons biopython seaborn r-coda
conda install -c conda-forge r-phangorn
conda activate phyload
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

## [`indices/`](indices)
Scripts needed to calculate test statistics used in the study.

### [`align_mi.py`](indices/align_mi.py)
Calculates all mutual information for _all_ pairs of sites in the alignment and prints to a file.
These are later summarized with <Insert method here> to turn into a single test statistic, like skewness.

To process a single alignment, use
```bash
$ <Insert cmd here>
```

### [`goldman_yang_1993.R`](indices/goldman_yang_1993.R)
Calculates the "unconstrained likelihood" test statistic of Goldman and Yang (1993).
This statistic is max(lnL(alignment)) under a multinomial model where all columns are iid from a multinomial distribution on the 4^n_taxa site patterns possible (assuming DNA).

Prints a file where the first line is "GY93" and the second line is the value of the test statistic for this alignment.

To process a single alignment, use
```bash
$ Rscript indices/goldman_yang_1993.R <path_to_aln> <output_file>
```

### [`pairwise_variance.R`](indices/pairwise_variance.R)
Calculates our new pairwise variance test statistic.
This statistic is var(all_pairwise_hamming_distances), the variance of the hamming distance (number of differences) for all pairwise comparisons of taxa in the alignment.

Prints a file where the first line is "PV" and the second line is the value of the test statistic for this alignment.

To process a single alignment, use
```bash
$ Rscript indices/pairwise_variance.R <path_to_aln> <output_file>
```
### [`singleton_fixed_tree.R`](indices/singleton_fixed_tree.R)
Calculates the p-value against the null hypothesis that all sites are iid (under some CTMC model of character evolution) on a given phylogeny.
For details of this test, see the paper.

Prints a file where the first line is "ST" and the second line is the p-value for this alignment.

To process a single alignment, use
```bash
$ Rscript indices/singleton_fixed_tree.R <path_to_aln> <path_to_tree> <output_file>
```

## [`simulation_scripts/`](simulation_scripts)
Scripts needed to simulate a single cell in the simulation study.

### [`simulate_alns.Rev`](simulation_scripts/simulate_alns.Rev)

To perform a single simulation, use
```bash
$ rb simulation_scripts/simulate_alns.Rev --args <n_iid> <n_epi> <d> <seed> <outbase> <treepath> <config>
```
Separate iid and epistatic alignments will be written in nexus format to `<outbase>/iid_aln.nex` and `<outbase>/epi_aln.nex`.
NOTE: The epistatic alignment are not DNA but characters 0-F, where 0=AA, 1=AC, ..., F=TT.
NOTE: The above assumes `rb` is on `$PATH`, if it is not replace `rb` with `<path/to/rb>`

The arguments are
- `n_iid`: number of iid sites
- `n_epi`: number of epistatically paired sites (must be even)
- `d`: *d* parameter for Nasrallah-Huelsenbeck model
- `seed`: random seed passed to RevBayes
- `outbase`: path to write output alignments
- `treepath`: path to treefile
- `config`: path to file with simulation parameters for all other (non-d) substitution models (iid and epistatic)

### Utilities

- [`merge_alns.R`](simulation_scripts/merge_alns.R)
Merge an iid nexus alignment `<iid_aln>` and and epistatic nexus alignment `<epi_aln>` (as produced by `mk_alns.R`) and write to a new nexus alignment `<merged_aln>` with
```bash
$ Rscript simulation_scripts/merge_alns.R <iid_aln> <epi_aln> <merged_aln>
```
This script will replace the 0-F characters in `<epi_aln>` with paired sites in `<merged_aln>`.
- [`epistatic_doublet_model_stub.Rev`](simulation_scripts/epistatic_doublet_model_stub.Rev) Core Rev script for the epistatic model.
- [`rev_model_template.Rev`](simulation_scripts/rev_model_template.Rev) A template for a Rev script to simulate an alignment in two parts, part epistatic and part purely iid sites, missing key parameters.
- [`simulation_tree.tre`](simulation_scripts/simulation_tree.tre) The treefile used for simulating alignments.
Rev points to this automatically for simulating.

## [`analysis_scripts/`](analysis_scripts)
Scripts needed to run a RevBayes analysis on a single cell in the simulation study.

### [`run_analysis.Rev`](analysis_scripts/run_analysis.Rev)

To analyze a single simulation, use
```bash
$ rb simulation_scripts/run_analysis.Rev --args <seed> <target_aln> <outbase>
```
Runs RevBayes on target_aln, stores output in format needed for posterior predictive simulation into `<outbase>/stochastic_variables.log` and `<outbase>/stochastic_variables_run_<1|2>.log`.

The arguments are
- `seed`: random seed passed to RevBayes
- `target_aln`: path to alignment to analyze
- `outbase`: path to write output files

### [`tree_distances.R`](analysis_scripts/tree_distances.R)

To compute a variety of summaries of distance from true tree to posterior trees, use
```bash
$ Rscript simulation_scripts/tree_distances.R <run1> <run2> <true_tree>
```
Computes the posterior distribution of RF and KF distances (that is, all RF and KF distances from trees sampled by the MCMC to the true tree).
Summarizes by taking the mean, median, min, max, and 2.5%, 5%, 95%, and 97.5% of the distribution.
Prints to stdout a tsv of the names of the summaries (row 1) and the summaries (row2).

The arguments are
- `run1`: path to posterior stochastic variables log for replicate/chain 1
- `run2`: path to posterior stochastic variables log for replicate/chain 2
- `true_tree`: path to the true tree

### [`split_based_metrics.R`](analysis_scripts/split_based_metrics.R)

To compute summaries related to the posterior distribution on splits, use
```bash
$ Rscript simulation_scripts/split_based_metrics.R <run1> <run2> <true_tree>
```
Computes how resolved the majority-rule consensus tree is (as a fraction of a fully bifurcating tree) and what percent of splits in the MRC tree are not in the true tree.
Prints to stdout a tsv of the names of the summaries (row 1) and the summaries (row2).

The arguments are
- `run1`: path to posterior stochastic variables log for replicate/chain 1
- `run2`: path to posterior stochastic variables log for replicate/chain 2
- `true_tree`: path to the true tree

### [`branch_length_measures.R`](analysis_scripts/branch_length_measures.R)

To compute a variety of summaries of tree and branch lengths, use
```bash
$ Rscript simulation_scripts/branch_length_measures.R <run1> <run2> <true_tree>
```
Computes the posterior distribution of tree length (sum of all branch lengths), sum of all tip branch lengths, and tree span (longest tip-to-tip distance).
Reports this relative to the true length, as a proportion (*e.g.* estimated_tree_length/true_tree_length).
Summarizes by taking the mean, median, min, max, and 2.5%, 5%, 95%, and 97.5% of the distribution.
Prints to stdout a tsv of the names of the summaries (row 1) and the summaries (row2).

The arguments are
- `run1`: path to posterior stochastic variables log for replicate/chain 1
- `run2`: path to posterior stochastic variables log for replicate/chain 2
- `true_tree`: path to the true tree

### [`run_PPS.R`](analysis_scripts/run_PPS.R)

To perform posterior predictive simulation based on the posterior from a single analysis, use
```bash
$ Rscript simulation_scripts/run_PPS.R <seed> <target_aln> <outbase> <rev_path>
```
Simulates alignments under the posterior predictive distribution into `<outbase>/PPS`.
Each alignment will appear in its own directory, `<outbase>/PPS/posterior_predictive_sim_i`.
Current setup will produce 102 simulated alignments per MCMC run.
We use an R wrapper around the script `run_pps.Rev` because we run PPS separately on each logfile and combine them, requiring some temporary directories lest we overwrite files.

The arguments are
- `seed`: random seed passed to RevBayes
- `target_aln`: path to alignment to analyze
- `outbase`: path to write output files
- `rev_path`: path to RevBayes executable

### Utilities

- [`diagnose_convergence.R`](analysis_scripts/diagnose_convergence.R)
Prints to stdout two convergence diagnostics for to filter out any analyses where MCMC convergence is suspect. Prints in order ASDSF (average standard deviation of split frequencies between two chains, split in half, ignoring splits of total frequency < 5%), PSRF (potential scale reduction factor between two chains, split in half, on the tree length). Run with,
```bash
$ Rscript analysis_scripts/diagnose_convergence.R <dir>
```
- [`analysis_template.Rev`](analysis_scripts/epistatic_doublet_model_stub.Rev) A template for a Rev script to analyze an alignment, missing key parameters.
