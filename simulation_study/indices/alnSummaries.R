library(phangorn)

# This is designed to be called with Rscript with several arguments
# arg1: path to alignment

# Call this script from phyload/simulation_study

# This writes to file a tsv of p-values

# Get arguments
args = commandArgs(trailingOnly=TRUE)

if (!length(args) == 1) {
  stop("This script requires 1 arguments")
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

out.dir <- dirname(args[1])

#' Compute the number of alignment sites in each of Lewis' 15 bins
#'
#' @param aln              phyDat alignment.
#' @param drop.ambiguous   TRUE/FALSE should columns with ambiguous characters be dropped?
#'
#' @details This function calculates the value of Goldman's (1993) "unconstrained likelihood" to use as a test statistic, as in Bollback (2002).
#' Specifically, the value computed is given by equation 6b in Goldman (1993) and equation 7 in Bollback (2002).
#'
#' @return The value of the test statistic.
countBins <- function(aln) {
  # recover()
  
  # For ease of grepping
  aln_mat <- tolower(as.character(aln))
  
  NCOL <- ncol(aln_mat)
  
  # Begin by getting list of what columns contain which of the four bases
  has.a.plus <- unique(which(aln_mat == "a",arr.ind=T)[,2])
  has.c.plus <- unique(which(aln_mat == "c",arr.ind=T)[,2])
  has.g.plus <- unique(which(aln_mat == "g",arr.ind=T)[,2])
  has.t.plus <- unique(which(aln_mat == "t",arr.ind=T)[,2])
  
  # Define columns that are missing each of the four bases
  not.a <- c(1:NCOL)[-has.a.plus]
  not.c <- c(1:NCOL)[-has.c.plus]
  not.g <- c(1:NCOL)[-has.g.plus]
  not.t <- c(1:NCOL)[-has.t.plus]
  
  # Check for gap-only columns
  has.none <- intersect(intersect(not.a,not.c),intersect(not.g,not.t))
  
  if ( length(has.none) > 0 ) {
    warning(paste0("There are ",length(has.none)," gap only columns in this alignment"))
    
    # # aln_mat <- aln_mat[,-has.none]
    # 
    # NCOL <- ncol(aln_mat)
    # 
    # # Begin by getting list of what columns contain which of the four bases
    # has.a.plus <- unique(which(aln_mat == "a",arr.ind=T)[,2])
    # has.c.plus <- unique(which(aln_mat == "c",arr.ind=T)[,2])
    # has.g.plus <- unique(which(aln_mat == "g",arr.ind=T)[,2])
    # has.t.plus <- unique(which(aln_mat == "t",arr.ind=T)[,2])
    # 
    # # Define columns that are missing each of the four bases
    # not.a <- c(1:NCOL)[-has.a.plus]
    # not.c <- c(1:NCOL)[-has.c.plus]
    # not.g <- c(1:NCOL)[-has.g.plus]
    # not.t <- c(1:NCOL)[-has.t.plus]
    
  }
  
  # Find columns with at least two of the bases
  has.a.c.plus <- intersect(has.a.plus,has.c.plus)
  has.a.g.plus <- intersect(has.a.plus,has.g.plus)
  has.a.t.plus <- intersect(has.a.plus,has.t.plus)
  has.c.g.plus <- intersect(has.c.plus,has.g.plus)
  has.c.t.plus <- intersect(has.c.plus,has.t.plus)
  has.g.t.plus <- intersect(has.g.plus,has.t.plus)
  
  # Find columns with all 4 bases
  has.a.c.g.t <- intersect(intersect(has.a.plus,has.c.plus),intersect(has.g.plus,has.t.plus))
  
  # Using sets, find columns with exactly three bases
  has.a.c.g <- intersect(intersect(has.a.c.plus,has.g.plus),not.t)
  has.a.c.t <- intersect(intersect(has.a.c.plus,has.t.plus),not.g)
  has.a.g.t <- intersect(intersect(has.a.g.plus,has.t.plus),not.c)
  has.c.g.t <- intersect(intersect(has.c.g.plus,has.t.plus),not.a)
  
  # Using sets, find columns with exactly two bases
  has.a.c <- intersect(intersect(has.a.c.plus,not.g),not.t)
  has.a.g <- intersect(intersect(has.a.g.plus,not.c),not.t)
  has.a.t <- intersect(intersect(has.a.t.plus,not.c),not.g)
  has.c.g <- intersect(intersect(has.c.g.plus,not.a),not.t)
  has.c.t <- intersect(intersect(has.c.t.plus,not.a),not.g)
  has.g.t <- intersect(intersect(has.g.t.plus,not.a),not.c)
  
  # Using everything we know, find columns with one base
  has.a <- intersect(intersect(not.c,not.g),not.t)
  has.c <- intersect(intersect(not.a,not.g),not.t)
  has.g <- intersect(intersect(not.a,not.c),not.t)
  has.t <- intersect(intersect(not.a,not.c),not.g)
  
  per.site <- numeric(NCOL)
  per.site[has.a] <- 1
  per.site[has.c] <- 2
  per.site[has.g] <- 3
  per.site[has.t] <- 4
  per.site[has.a.c] <- 5
  per.site[has.a.g] <- 6
  per.site[has.a.t] <- 7
  per.site[has.c.g] <- 8
  per.site[has.c.t] <- 9
  per.site[has.g.t] <- 10
  per.site[has.a.c.g] <- 11
  per.site[has.a.c.t] <- 12
  per.site[has.a.g.t] <- 13
  per.site[has.c.g.t] <- 14
  per.site[has.a.c.g.t] <- 15
  per.site[has.none] <- -Inf
  
  return(per.site)
}

site.cats <- countBins(aln)
all.bins <- sapply(1:15,function(i){sum(site.cats == i)})
cat("LB","\n",paste0(all.bins,collapse=","),sep="",file=paste0(out.dir,"/LB.txt"))

prop.inv <- sum(all.bins[1:4])/sum(all.bins)
cat("PI","\n",prop.inv,sep="",file=paste0(out.dir,"/PI.txt"))

site.nucs <- sapply(site.cats,function(s){
  if (s %in% c(1:4)) {
    return(1)
  } else if (s %in% c(5:10)) {
    return(2)
  } else if (s %in% c(11:14)) {
    return(3)
  } else if (s == 15) {
    return(4)
  } else {
    return(NA)
  }
})

bcd <- mean(site.nucs)
cat("BCD","\n",bcd,sep="",file=paste0(out.dir,"/BCD.txt"))

bcv <- var(site.nucs)
cat("BCV","\n",bcv,sep="",file=paste0(out.dir,"/BCV.txt"))
