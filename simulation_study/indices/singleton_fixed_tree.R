library(phangorn)

# This is designed to be called with Rscript with several arguments
# arg1: path to alignment
# arg2: path to treefile (single tree or posterior)
# arg3: path to desired output file

# Call this script from phyload/simulation_study

# This writes to file a tsv of p-values

# Get arguments
args = commandArgs(trailingOnly=TRUE)

if (!length(args) == 3) {
  stop("This script requires 3 arguments")
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
phy <- read.tree(args[2])

out.file <- args[3]

# Singleton test for independence of sites.
#'
#' @param aln              phyDat alignment.
#' @param phy              phylo-style tree or multiPhylo-style list of trees
#' @param nrep             the test uses Monte Carlo simulation with nrep simulations to calculate the p-value, default 2000 as in R's chisq.test
#' @param min.edge.length  if any branch lengths are below this value (that is, they are essentially 0 either due to rounding or inference programs thresholding branches to 0) we set the branch to this length to avoid dividing by 0
#' 
#' @details This function tests if sites are independent using singleton sites only.
#'          If there is a posterior distribution on trees, this will perform the test on all trees
#' 
#' @return The valuep-value against independence.

testSiteIID <- function(aln,phy,nrep=2000,min.edge.length=1e-8) {
  # recover()
  
  # If phy is a single tree, make phy a multiPhylo so we can use loop regardless
  if ( class(phy) == "phylo" ) {
    tmp <- list(phy)
    class(tmp) <- "multiPhylo"
    phy <- tmp
  } else if ( class(phy) == "list" ) {
    if ( class(phy[[1]]) == "phylo" ) {
      class(phy) <- "multiPhylo"
    } else {
      return(NA)
    }
  } else if (class(phy) != "multiPhylo") {
    return(NA)
  }
  
  aln_num <- do.call(rbind,aln)
  
  ordered_taxa <- sort(row.names(aln_num))
  
  aln_tab <- apply(aln_num,2,table)
  
  ntaxa <- length(aln)
  
  site_pattern_counts <- attributes(aln)$weight
  
  # Reduce alignment to count data
  # singletons is a character vector containing one occurrence each time that taxon is the singleton in a singleton site
  possible_singletons <- which(lengths(aln_tab) == 2)
  singletons <- lapply(possible_singletons,function(i){
    x <- aln_tab[[i]]
    # A singleton is a site with two kinds of bases with one present only once
    if ( any(x == 1) ) {
      the_base <- as.numeric(names(x[x == 1]))
      return(rep(row.names(aln_num)[which(aln_num[,i] == the_base)],site_pattern_counts[i]))
    } else {
      return(NULL)
    }
  })
  # Make singletons an integer vector by using taxon number instead of name
  singletons <- unlist(singletons)
  singletons <- sapply(singletons,function(x){which(ordered_taxa == x)})

  # Table of the number of observed singletons for each taxon
  obs <- sapply(1:ntaxa,function(i){sum(singletons == i)})
  
  if (sum(obs) == 0) {
    stop("No singleton sites in this alignment")
  }
  
  # p-values against independence
  p <- numeric(length(phy))
  
  for (n in 1:length(phy)) {
    tree <- phy[[n]]
    tree_tip_edge_lengths <- sapply(1:ntaxa,function(i){
      tip_index <- which(tree$tip.label == ordered_taxa[i])
      tree$edge.length[which(tree$edge[,2] == tip_index)]
    })
    tree_tip_edge_lengths[tree_tip_edge_lengths < min.edge.length] <- min.edge.length
    
    tree_singleton_probs <- tree_tip_edge_lengths/sum(tree_tip_edge_lengths)
    
    # Use Monte Carlo simulation to calculate p-value
    # p[n] <- chisq.test(obs,p=tree_singleton_probs,simulate.p.value = T)$p.value
    
    # This is faster than using R's built in version, but equivalent
    n_singletons <- length(singletons)
    expected <- n_singletons * tree_singleton_probs
    xi_sq <- sum( ((obs - expected)^2)/expected )
    
    sim <- rmultinom(nrep,n_singletons,tree_singleton_probs)
    sim <- (sim - expected)^2/expected
    null_dist <- colSums(sim)
    
    # (n_greater + 1) / (n + 1) is how R calculates the p-value
    p[n] <- (sum(null_dist > xi_sq)+1)/(nrep+1)
  }
  return(p)
}

p.val <- testSiteIID(aln,phy)

cat("ST","\n",p.val,sep="",file=out.file)
