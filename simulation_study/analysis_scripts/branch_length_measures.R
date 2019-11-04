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

all.trees <- c(run1,run2)

ntax <- length(all.trees[[1]]$tip.label)

# tree length (sum of all branches)
all.tl <- unlist(lapply(all.trees,function(phy){sum(phy$edge.length)}))
# total length of tip branches
all.tl.tip <- unlist(lapply(all.trees,function(phy){sum(phy$edge.length[phy$edge[,2] <= ntax])}))
# longest tip-to-tip distance
all.span <- unlist(lapply(all.trees,function(phy){max(cophenetic.phylo(phy))}))

tl.min    <- min(all.tl)
tl.max    <- max(all.tl)
tl.mean   <- mean(all.tl)
tl.median <- median(all.tl)
tl.q025   <- quantile(all.tl,0.025)
tl.q05    <- quantile(all.tl,0.05)
tl.q95    <- quantile(all.tl,0.95)
tl.q975   <- quantile(all.tl,0.975)

tl.tip.min    <- min(all.tl.tip)
tl.tip.max    <- max(all.tl.tip)
tl.tip.mean   <- mean(all.tl.tip)
tl.tip.median <- median(all.tl.tip)
tl.tip.q025   <- quantile(all.tl.tip,0.025)
tl.tip.q05    <- quantile(all.tl.tip,0.05)
tl.tip.q95    <- quantile(all.tl.tip,0.95)
tl.tip.q975   <- quantile(all.tl.tip,0.975)

span.min    <- min(all.span)
span.max    <- max(all.span)
span.mean   <- mean(all.span)
span.median <- median(all.span)
span.q025   <- quantile(all.span,0.025)
span.q05    <- quantile(all.span,0.05)
span.q95    <- quantile(all.span,0.95)
span.q975   <- quantile(all.span,0.975)

tl.summaries <- c(tl.min,tl.max,tl.mean,tl.median,tl.q025,tl.q05,tl.q95,tl.q975)
tl.tip.summaries <- c(tl.tip.min,tl.tip.max,tl.tip.mean,tl.tip.median,tl.tip.q025,tl.tip.q05,tl.tip.q95,tl.tip.q975)
span.summaries <- c(span.min,span.max,span.mean,span.median,span.q025,span.q05,span.q95,span.q975)

tl.summary.names <- c("tl.min","tl.max","tl.mean","tl.median","tl.q025","tl.q05","tl.q95","tl.q975")
tl.tip.summary.names <- c("tl.tip.min","tl.tip.max","tl.tip.mean","tl.tip.median","tl.tip.q025","tl.tip.q05","tl.tip.q95","tl.tip.q975")
span.summary.names <- c("span.min","span.max","span.mean","span.median","span.q025","span.q05","span.q95","span.q975")

cat(c(tl.summary.names,tl.tip.summary.names,span.summary.names), file=stdout(), sep="\t")
cat("\n",sep="")
cat(c(tl.summaries,tl.tip.summaries,span.summaries), file=stdout(), sep="\t")
