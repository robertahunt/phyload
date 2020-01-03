library(phangorn)
library(coda)

# This is designed to be called with Rscript with several arguments
# arg1: path to posterior sample 1 .tree file
# arg2: path to posterior sample 2 .tree file

# Call this script from phyload/simulation_study

# This prints two convergence diagnostics (ASDSF and tree length PSRF) to stdout

# Get arguments
args = commandArgs(trailingOnly=TRUE)

if (!length(args) == 2) {
  stop("This script requires 2 arguments")
}

run1 <- read.tree(args[1])
run2 <- read.tree(args[2])


# Estimates the average standard deviation of split frequencies for multiple chains of trees
# Arguments:
#   chains: all chains as a list (need not be multiphylo class)
#   split: should we split each chain in half before calculation (adds within chain component to overall convergence)
#   min.freq.cutoff: we will ignore all splits that occur at less than this frequency overall
# Returns: the ASDSF
ASDSF <- function(chains,split.chains=FALSE,min.freq.cutoff=0.05) {
  # recover()

  all_trees <- do.call(c,chains)
  class(all_trees) <- "multiPhylo"

  # We want to loop over all trees, and for any tree know what chain it belongs to
  if (split.chains) {
    which_chain <- lapply(1:length(chains),function(i){
      n_1 <- floor(length(chains[[i]])/2)
      n_2 <- length(chains[[i]]) - n_1
      c(rep(2*(i-1)+1,n_1),rep(2*i,n_2))
    })
    which_chain <- unlist(which_chain)
  } else {
    which_chain <- lapply(1:length(chains),function(i){rep(i,length(chains[[i]]))})
    which_chain <- unlist(which_chain)
  }

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
  counts <- matrix(0,nrow=(ntax-3)*length(all_trees),ncol=length(chains)*ifelse(split.chains,2,1))

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
      counts[(n_seen+1):(n_seen+sum(!seen_these_already)),which_chain[i]] <- counts[(n_seen+1):(n_seen+sum(!seen_these_already)),which_chain[i]] + 1
    }

    # Add seen splits to count for appropriate chain
    counts[1:n_seen,which_chain[i]][saw_them_here] <- counts[1:n_seen,which_chain[i]][saw_them_here] + 1

    n_seen <- n_seen + sum(!seen_these_already)

    count_sum[i] <- sum(counts)
  }

  # Remove unneeded parts of counts
  counts <- counts[1:n_seen,]

  # Make counts into probabilities
  for (i in 1:(length(chains)*ifelse(split.chains,2,1))) {
    counts[,i] <- counts[,i]/sum(which_chain == i)
  }

  sdsf <- apply(counts,1,sd)

  above_cutoff <- apply(counts,1,function(p){mean(p) > min.freq.cutoff})

  return(mean(sdsf[above_cutoff]))

}

# Estimates the PSRF of the lengths of the trees
# Arguments:
#   chains: all chains as a list (need not be multiphylo class)
# Returns: the PSRF
treeLengthPSRF <- function(chains,split.chains=FALSE) {
  # Get all tree lengths, preserve chains
  chains <- lapply(chains,function(chain){
    unlist(lapply(chain,function(phy){sum(phy$edge.length)}))
  })

  n <- lengths(chains)
  m <- length(chains)

  # We're assuming these are all runs of equal length of the same model, all MCMC outputs should be the same size
  if ( length(unique(n)) != 1 ) {
    stop("Input chains are not of the same length")
  }
  n <- n[1]
  # recover()

  if (split.chains) {
    m <- 2*length(chains)

    if ( n %% 2 == 1 ) {
      # Chains have odd numbers of samples, n will be variable and we must assign the odd sample to one half of each chain
      # warning("Chains have odd numbers of samples, discarding last sample in calculation of diagnostics")
      chains <- lapply(chains,function(x){x[-length(x)]})
      n <- (n-1)/2
    } else {
      n <- n/2
    }

    # Split chains in half, make sure each half is a matrix
    short_chains <- vector("list",m)
    index <- 0
    for (i in 1:(m/2)) {
      index <- index + 1
      short_chains[[index]] <- as.mcmc(chains[[i]][1:(n)])
      index <- index + 1
      short_chains[[index]] <- as.mcmc(chains[[i]][(n+1):(2*n)])
    }
    chains <- short_chains
  }

  chains <- lapply(chains,as.mcmc)

  gelman.diag(chains)$psrf[1]

}

asdsf <- ASDSF(list(run1, run2))
psrf  <- treeLengthPSRF(list(run1, run2))

cat("asdsf", "psrf\n", file=stdout(), sep="\t")
cat(asdsf, psrf, file=stdout(), sep="\t")
