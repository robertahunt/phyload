library(phangorn)

# This is designed to be called with Rscript with several arguments
# arg1: seed for Rev
# arg2: alignment that PPS is to be run on
# arg3: path to write output files
# arg4: command for calling RevBayes

# Call this script from phyload/simulation_study

# Get arguments
args = commandArgs(trailingOnly=TRUE)

if (!length(args) == 4) {
  stop("This script requires 4 arguments")
}

seed   <- args[1]
aln    <- args[2]
outdir <- args[3]
rb     <- args[4]

pps.dir <- paste0(outdir,"/PPS")

# Run PPS on each logfile separately in temporary subdirectories so that we don't overwrite PPS alignments
dir1 <- paste0(outdir,1)
dir2 <- paste0(outdir,2)

dir.create(dir1)
dir.create(dir2)
dir.create(pps.dir)

system2("rb",args=c("analysis_scripts/run_pps.Rev","--args",seed,aln,outdir,1))
system2("rb",args=c("analysis_scripts/run_pps.Rev","--args",seed,aln,outdir,2))

n_pps <- length(list.files(paste0(dir1,"/PPS")))

# Merge all PPS and re-number second logfile runs to come after first log
system2("mv",args=c(paste0(dir1,"/PPS/*"),paste0(pps.dir,"/")))

for (i in 1:n_pps) {
  system2("mv",args=c(paste0(dir2,"/PPS/posterior_predictive_sim_",i),paste0(dir2,"/PPS/posterior_predictive_sim_",i+n_pps)))
}

system2("mv",args=c(paste0(dir2,"/PPS/*"),paste0(pps.dir,"/")))

# Remove temp dirs
file.remove(paste0(dir1,"/PPS"))
file.remove(dir1)

file.remove(paste0(dir2,"/PPS"))
file.remove(dir2)
