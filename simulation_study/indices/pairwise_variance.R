library(phangorn)

# This is designed to be called with Rscript with several arguments
# arg1: path to alignment (or folder with posterior predictive alignments)
# arg2: path to desired output file

# Call this script from phyload/simulation_study

# This writes to file a tsv of p-values

# Get arguments
args = commandArgs(trailingOnly=TRUE)

if (!length(args) == 2) {
  stop("This script requires 2 arguments")
}

file.format <- strsplit(basename(args[1]),".",fixed=TRUE)[[1]]
file.format <- tolower(file.format[length(file.format)])
if ( grepl("nex",file.format) ) {
  file.format <- "nexus"
} else if ( grepl("phy",file.format) ) {
  file.format <- "phylip"
} else if ( grepl("fa",file.format) ) {
  file.format <- "fasta"
} else {
  stop("Unrecognized file format of alignment.")
}
aln <- read.phyDat(args[1],format=file.format)

out.file <- args[2]

#' Compute the posterior predictive test statistic of Goldman and Yang (1993)
#'
#' @param aln              phyDat alignment.
#' 
#' @details This function calculates the variance of the number of substitutions in all pairwise comparisons.
#' 
#' @return The value of the test statistic.

pairwiseVariance <- function(aln) {
  pairs <- lapply(1:(length(aln)-1),function(i){
    lapply((i+1):length(aln),function(j){
      sum(aln[[i]] != aln[[j]])
    })
  })
  return(var(unlist(pairs)))
}

t.stat <- pairwiseVariance(aln)

cat("PV","\n",t.stat,sep="",file=out.file)
