library(phangorn)

# This is designed to be called with Rscript with several arguments
# arg1: path to posterior trees from run 1
# arg2: path to posterior trees from run 2
# arg3: path to true tree

# Call this script from phyload/simulation_study

# This prints several 

# Get arguments
args = commandArgs(trailingOnly=TRUE)

if (!length(args) == 3) {
  stop("This script requires 3 arguments")
}

run1 <- read.tree(args[1])
run2 <- read.tree(args[2])

true <- read.tree(args[3])


# Finds all resolved (bifurcating/non-polytomy) nodes in the MRC tree
# Arguments:
#   chains: all chains as a list (need not be multiphylo class)
# Returns: the resolved splits in the MRC tree
findResolvedSplitsMRC <- function(chains) {
  # recover()
  
  all_trees <- do.call(c,chains)
  class(all_trees) <- "multiPhylo"
  
  taxa <- all_trees[[1]]$tip.label
  
  ntax <- length(taxa)
  
  # To avoid constantly making new rows, we upper bound the number of unique splits
  # An unrooted tree has 2ntax-3 internal branches, each of which defines a split
  # ntax-3 of these are non-trivial splits
  # Thus there are no more than (ntax-3)*length(all_trees) unique splits
  all_splits <- character((ntax-3)*length(all_trees))
  
  # We also track the number of seen splits, to avoid searching in NAs
  n_seen <- 0
  
  # To count occurences in each chain
  counts <- numeric((ntax-3)*length(all_trees))
  
  count_sum <- numeric((ntax-3)*length(all_trees))
  
  # Store splits alphabetically
  # colnames(all_splits) <- sort(taxa)
  
  for (i in 1:length(all_trees)) {
    # splits objects are annoying
    these_splits <- as.matrix(as.splits(all_trees[[i]]))
    
    # alphabetize
    these_splits <- these_splits[,order(colnames(these_splits))]
    
    # remove trivial splits (only one taxon or all taxa)
    trivial <- rowSums(these_splits) == 1 | rowSums(these_splits) == ntax
    
    these_splits <- these_splits[!trivial,]
    
    # Wooo what could be more fun than handling strings in R?
    these_splits <- apply(these_splits,1,paste0,collapse="")
    # seen_these_already <- these_splits %in% all_splits[1:n_seen]
    # saw_them_here <- all_splits[1:n_seen] %in% these_splits
    
    # We need to check if we've seen either the split or its complement
    repolarized <- these_splits
    repolarized <- gsub("1","x",repolarized)
    repolarized <- gsub("0","1",repolarized)
    repolarized <- gsub("x","0",repolarized)
    
    # seen_these_already[!seen_these_already] <- repolarized[!seen_these_already] %in% all_splits[1:n_seen]
    # saw_them_here[!seen_these_already]
    seen_these_already <- (these_splits %in% all_splits[1:n_seen]) | (repolarized %in% all_splits[1:n_seen])
    saw_them_here <- (all_splits[1:n_seen] %in% these_splits) | (all_splits[1:n_seen] %in% repolarized)
    
    # Add new splits to master list and count them, if there are any new splits
    if (sum(!seen_these_already) > 0) {
      all_splits[(n_seen+1):(n_seen+sum(!seen_these_already))] <- these_splits[!seen_these_already]
      counts[(n_seen+1):(n_seen+sum(!seen_these_already))] <- counts[(n_seen+1):(n_seen+sum(!seen_these_already))] + 1
    }
    
    # Add seen splits to count for appropriate chain
    counts[1:n_seen][saw_them_here] <- counts[1:n_seen][saw_them_here] + 1
    
    n_seen <- n_seen + sum(!seen_these_already)
  }
  
  # Remove unneeded parts
  counts <- counts[1:n_seen]
  all_splits <- all_splits[1:n_seen]
  
  # Make counts into probabilities
  counts <- counts/length(all_trees)
  
  return(all_splits[counts > 0.5])
  
}

# Finds all splits in a list that are not in the true tree
# Arguments:
#   true.tree: the true tree
#   splits: a reference list of splits, as characters (e.g. "100111000")
# Returns: the number of resolved nodes in the MRC tree
wrongSplits <- function(true.tree,splits) {
  # recover()
  
  ntax <- length(true.tree$tip.label)
  
  # Get splits in true tree
  true_splits <- as.matrix(as.splits(true))
  true_splits <- true_splits[,order(colnames(true_splits))]
  trivial <- rowSums(true_splits) == 1 | rowSums(true_splits) == ntax
  true_splits <- true_splits[!trivial,]
  true_splits <- apply(true_splits,1,paste0,collapse="")
  
  # correct splits given current polarization
  correct_as_is <- splits %in% true_splits
  
  # We need to check if we've seen either the split or its complement
  repolarized <- splits
  repolarized <- gsub("1","x",repolarized)
  repolarized <- gsub("0","1",repolarized)
  repolarized <- gsub("x","0",repolarized)
  
  correct_repolarized <- repolarized %in% true_splits
  
  correct <- correct_as_is | correct_repolarized
  
  return(splits[!correct])
}

ntaxa <- length(true$tip.label)

# Get all splits in MRC from posterior and calculate what percent resolved it is
resolved <- findResolvedSplitsMRC(list(run1,run2))

percent_resolved <- length(resolved)/(ntaxa-3)

resolved_wrong <- wrongSplits(true,resolved)

percent_resolved_wrong <- length(resolved_wrong)/length(resolved)

cat("mrc_percent_resolved", "mrc_percent_wrong_splits",file=stdout(), sep="\t")
cat("\n")
cat(percent_resolved, percent_resolved_wrong, file=stdout(), sep="\t")
