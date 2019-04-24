# Function to pool alignments
source("simulation_scripts/merge_alignments.R")

# This is designed to be called with Rscript with several arguments
# arg1: the (relative) file path to directory we're simulating in
# arg2: do we want to delete the component alignments?

# Call this script from phyload/simulation_study

# Get arguments
args = commandArgs(trailingOnly=TRUE)

if (!length(args) == 2) {
  stop("This script requires 2 arguments")
}

out.dir    <- args[1]
delete     <- args[2]

replicate <- basename(dirname(out.dir))

iid.aln  <- paste0(out.dir,"/base_alignment_",replicate,".nex")
epi.aln  <- paste0(out.dir,"/epistatic_alignment_",replicate,".nex")
out.file <- paste0(out.dir,"/",replicate,".nex")

# Call merge function
mergeAlignments(iid.aln,epi.aln,out.file)

if (delete) {
  system2("rm",iid.aln)
  system2("rm",epi.aln)
}