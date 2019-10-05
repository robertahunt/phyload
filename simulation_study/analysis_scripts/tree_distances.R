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

post <- c(run1,run2)
class(post) <- "multiPhylo"

all.rf <- RF.dist(post,true)
all.kf <- RF.dist(post,true)

rf.min    <- min(all.rf)
rf.max    <- max(all.rf)
rf.mean   <- mean(all.rf)
rf.median <- median(all.rf)
rf.q025   <- quantile(all.rf,0.025)
rf.q05    <- quantile(all.rf,0.05)
rf.q95    <- quantile(all.rf,0.95)
rf.q975   <- quantile(all.rf,0.975)

kf.min    <- min(all.kf)
kf.max    <- max(all.kf)
kf.mean   <- mean(all.kf)
kf.median <- median(all.kf)
kf.q025   <- quantile(all.kf,0.025)
kf.q05    <- quantile(all.kf,0.05)
kf.q95    <- quantile(all.kf,0.95)
kf.q975   <- quantile(all.kf,0.975)

rf.summaries <- c(rf.min,rf.max,rf.mean,rf.median,rf.q025,rf.q05,rf.q95,rf.q975)
kf.summaries <- c(kf.min,kf.max,kf.mean,kf.median,kf.q025,kf.q05,kf.q95,kf.q975)

rf.summary.names <- c("rf.min","rf.max","rf.mean","rf.median","rf.q025","rf.q05","rf.q95","rf.q975")
kf.summary.names <- c("kf.min","kf.max","kf.mean","kf.median","kf.q025","kf.q05","kf.q95","kf.q975")

cat(c(rf.summary.names,kf.summary.names), file=stdout(), sep="\t")
cat("\n",sep="")
cat(c(rf.summaries,kf.summaries), file=stdout(), sep="\t")
