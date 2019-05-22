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

aln <- read.phyDat(args[1])
phy <- read.tree(args[2])

out.file <- args[3]

# Singleton test for independence of sites.
#'
#' @param aln              phyDat alignment.
#' @param phy              phylo-style tree or multiPhylo-style list of trees
#' 
#' @details This function tests if sites are independent using singleton sites only.
#'          If there is a posterior distribution on trees, this will perform the test on all trees
#' 
#' @return The valuep-value against independence.

testSiteIID <- function(aln,phy) {
  
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
  
  aln_char <- tolower(as.character(aln))
  
  ordered_taxa <- sort(row.names(aln_char))
  
  aln_tab <- apply(aln_char,2,table)
  
  ntaxa <- length(aln)
  
  # Reduce alignment to count data
  # singletons is a character vector containing one occurrence each time that taxon is the singleton in a singleton site
  singletons <- lapply(1:dim(aln_char)[2],function(i){
    x <- aln_tab[[i]]
    # A singleton is a site with two kinds of bases with one present only once
    if (length(x) == 2 && any(x) == 1) {
      the_base <- names(x[x == 1])
      return(row.names(aln_char)[which(aln_char[,i] == the_base)])
    } else {
      return(NULL)
    }
  })
  # Make singletons an integer vector by using taxon number instead of name
  singletons <- unlist(singletons)
  singletons <- sapply(singletons,function(x){which(ordered_taxa == x)})
  names(singletons) <- NULL
  
  # Table of the number of observed singletons for each taxon
  obs <- sapply(1:ntaxa,function(i){sum(singletons == i)})
  
  # p-values against independence
  p <- numeric(length(phy))
  
  for (n in 1:length(phy)) {
    tree <- phy[[n]]
    tree_tip_edge_lengths <- sapply(1:ntaxa,function(i){
      tip_index <- which(tree$tip.label == ordered_taxa[i])
      tree$edge.length[which(tree$edge[,2] == tip_index)]
    })
    
    tree_singleton_probs <- tree_tip_edge_lengths/sum(tree_tip_edge_lengths)
    
    expected <- length(singletons) * tree_singleton_probs
    xi_sq <- sum( ((obs - expected)^2)/expected )
    p[n] <- pchisq(xi_sq,ntaxa-1,lower.tail=FALSE)
  }
  return(p)
}

p.vals <- testSiteIID(aln,phy)

cat(p.vals,sep="\t",file=out.file)
