library(phangorn)

# This is designed to be called with Rscript with several arguments
# arg1: path to alignment
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
#' @param drop.ambiguous   TRUE/FALSE should columns with ambiguous characters be dropped?
#'
#' @details This function calculates the value of Goldman's (1993) "unconstrained likelihood" to use as a test statistic, as in Bollback (2002).
#' Specifically, the value computed is given by equation 6b in Goldman (1993) and equation 7 in Bollback (2002).
#'
#' @return The value of the test statistic.

unconstrainedLikelihood <- function(aln,drop.ambiguous=FALSE) {
  # Nomenclature follows Bollback's (2002) description of the test statistic.

  # For ease of grepping
  aln_mat <- tolower(as.character(aln))

  # Collapse sites into character strings, use table to compute number of occurences
  N_Xi <- table(apply(aln_mat,2,paste,collapse=""))

  if ( drop.ambiguous ) {
    has_ambiguity <- ceiling(which(!(N_Xi %in% c("a","c","g","t")))/length(N_Xi))
    N_Xi <- N_Xi[!has_ambiguity]
  }

  # N is the total number of sites
  N <- sum(N_Xi)

  # Compute statistic
  T_of_X <- sum(N_Xi * log(N_Xi)) - N*log(N)
  return(T_of_X)
}

t.stat <- unconstrainedLikelihood(aln)

cat("G93","\n",t.stat,sep="",file=out.file)
