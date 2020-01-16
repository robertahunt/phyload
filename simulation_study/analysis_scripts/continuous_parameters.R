# This is designed to be called with Rscript with several arguments
# arg1: path to posterior sample 1 stochastic variables trace
# arg2: path to posterior sample 2 stochastic variables trace

# Call this script from phyload/simulation_study

# This prints two convergence diagnostics (ASDSF and tree length PSRF) to stdout

# Get arguments
args = commandArgs(trailingOnly=TRUE)

if (!length(args) == 2) {
  stop("This script requires 2 arguments")
}

run1 <- read.table(args[1],header=TRUE,stringsAsFactors=FALSE)
run2 <- read.table(args[1],header=TRUE,stringsAsFactors=FALSE)

alpha <- c(run1$alpha,run2$alpha)

er <- c(run1$er,run2$er)
er <- gsub("[","",er,fixed=TRUE)
er <- gsub("]","",er,fixed=TRUE)
er <- lapply(er,function(x){
  as.numeric(strsplit(x,",")[[1]])
})
er <- do.call(rbind,er)

bf <- c(run1$pi,run2$pi)
bf <- gsub("[","",bf,fixed=TRUE)
bf <- gsub("]","",bf,fixed=TRUE)
bf <- lapply(bf,function(x){
  as.numeric(strsplit(x,",")[[1]])
})
bf <- do.call(rbind,bf)

all.par <- cbind(alpha,er,bf)
colnames(all.par) <- c("ASRF_gamma_shape","rate[AC]","rate[AG]","rate[AT]","rate[CG]","rate[CT]","rate[GT]","pi[A]","pi[C]","pi[G]","pi[T]")

cat(paste0(colnames(all.par),sep="\t"),"\n", file=stdout(), sep="")
cat(paste0(colMeans(all.par),sep="\t"),"\n", file=stdout(), sep="")
