#' Compute the posterior predictive test statistic of Bollback (2002).
#'
#' @param aln.in     Character. An alignment in matrix form, as read in with readAlignment().
#' 
#' @details This function calculates the value of Goldman's (1993) "unconstrained likelihood" to use as a test statistic, as in Bollback (2002).
#' Specifically, the value computed is given by equation 6b in Goldman (1993) and equation 7 in Bollback (2002).
#' Currently, all columns in the alignment containing gaps are discarded in this calculation.
#' 
#' @references 
#' Goldman, N. (1993). Statistical Tests of Models of DNA Substitution. Journal of Molecular Evolution, 36, 182-198.
#' 
#' Bollback, J. P. (2002). Bayesian Model Adequacy and Choice in Phylogenetics. Molecular Biology and Evolution, 19, 1171-1182.
#' 
#' @return The value of the test statistic.

unconstrainedLikelihood <- function(aln.mat) {
  # Nomenclature follows Bollback's (2002) description of the test statistic.
  
  # For ease of grepping
  aln_mat <- apply(aln.mat,1,tolower)
  
  # Collapse sites into character strings, use table to compute number of occurences
  N_Xi <- table(apply(aln_mat,2,paste,collapse=""))
  
  # Drop columns with gaps and missing data
  if (any(grepl("-",N_Xi,fixed=TRUE))) {
    N_Xi <- N_Xi[-which(grep("-",N_Xi,fixed=TRUE))]
  }
  if (any(grepl("?",N_Xi,fixed=TRUE))) {
    N_Xi <- N_Xi[-which(grep("?",N_Xi,fixed=TRUE))]
  }
  
  # N is the total number of sites
  N <- sum(N_Xi)
  
  # Compute statistic
  T_of_X <- sum(N_Xi * log(N_Xi)) - N*log(N)
  return(T_of_X)
}


#' Compute the posterior predictive test statistic of Huelsenbeck et al. (2001).
#'
#' @param aln.in     Character. An alignment in matrix form, as read in with readAlignment().
#' 
#' @details This function calculates the value of Goldman's (1993) posterior predictive test statistic for an alignment.
#' Specifically, the value computed is given by equation 6b in Goldman (1993) and equation 7 in Bollback (2002).
#' Currently, all sites with gaps are discarded in this calculation.
#' 
#' @references 
#' Huelsenbeck, J. P., Ronquist, F., Nielsen, R., and Bollback, J. P. (2001). Bayesian Inference of Phylogeny and Its Impact on Evolutionary Biology. Science, 294, 2310-2314.
#' 
#' @return The value of the test statistic.

chiSquaredHomogeneity <- function(aln.mat) {
  # recover()
  # For ease of grepping
  aln_mat <- apply(aln.mat,1,tolower)
  
  # Calculate useful information
  ntax <- nrow(aln.mat)
  bf_counts <- getCountsDistributeUncertainties(as.character(aln_mat))
  avg_base_freqs <- bf_counts/sum(bf_counts)
  
  species_deviances <- numeric(ntax)
  for (i in 1:ntax) {
    counts <- getCountsDistributeUncertainties(aln.mat[i,])
    species_base_freqs <- counts / sum(counts)
    species_deviances[i] <- sum( ((species_base_freqs-avg_base_freqs)^2) / avg_base_freqs )
  }
  
  # Compute statistic
  chi_squared <- sum(species_deviances)
  return(chi_squared)
}

#' Compute the mean biochemical diversity posterior predictive test statistic of Lartillot et al. (2007).
#'
#' @param aln.in     Character. An alignment in matrix form, as read in with readAlignment().
#' 
#' @details This function calculates the value of Lartillot et al.'s (2007) "mean biochemical diversity" to use as a test statistic.
#' Specifically, the value is the mean number of distinct nucleotides (A, C, G, and T) at each site.
#' Currently, gaps, missing data, and ambiguities are ignored.
#' 
#' @references 
#' Lartillot, N., Brinkmann, H., and Philippe, H. (2007). Suppression of long-branch attraction artefacts in the animal phylogeny using a site-heterogeneous model. BMC Evolutionary Biology, 7, 1-14.
#' 
#' @return The value of the test statistic.

biochemicalDiversity <- function(aln.mat) {
  # recover()
  # For ease of grepping, we make the alignment columns into lowercase strings
  aln_mat <- tolower(aln.mat)
  aln_mat <- apply(aln_mat,2,paste,collapse="")
  
  # calculate number of distinct nucleotides per site
  sitewise_diversity <- sapply(aln_mat,function(column) {
    # is each of A, C, G, and T present?
    has.base <- c(grepl("a",column),grepl("c",column),grepl("g",column),grepl("t",column))
    return(length(which(has.base)))
  })
  return(mean(sitewise_diversity))
}

#' Compute the variance of the site-soecific empirical frequencies of the bases.
#'
#' @param aln.in     Character. An alignment in matrix form, as read in with readAlignment().
#' 
#' @details This function calculates the value of Lartillot et al.'s (2007) "mean biochemical diversity" to use as a test statistic.
#' Specifically, the value is the mean number of distinct nucleotides (A, C, G, and T) at each site.
#' Currently, gaps, missing data, and ambiguities are ignored.
#' 
#' @references 
#' Feuda et al. (2017). Improved modeling of compositional heterogeneity supports sponges as sister to all other animals. Cell, 27, 3864-3870.
#' 
#' @return The value of the test statistic.

CATVar <- function(aln.mat) {
  # recover()
  # For ease of grepping, we make the alignment columns into lowercase strings
  aln_mat <- apply(aln.mat,1,tolower)
  
  # calculate number of distinct nucleotides per site
  base_counts <- apply(aln_mat,2,getCountsDistributeUncertainties)
  base_freqs <- base_counts/colSums(base_counts)
  
  # Variances
  base_freq_vars <- apply(base_freqs,1,var)
  
  return(mean(base_freq_vars))
}


#' Calculates y. as described in Lewis (2014).
#'  
#' @param aln.mat    Character. An alignment in matrix form, as read in with readAlignment().
#' 
#' @details This function calculates y. as in Lewis (2014), where y. is a vector of pattern counts.
#' This particular function defines either 7 or 15 site patterns.
#' The 15 site patterns are:
#' 
#'     1-4) Columns containing only one base (A, C, G, or T), not counting gaps, missing data, or other ambiguities
#'     
#'     5-10) Columns containing exactly two specific bases (AG, CT, etc, again, not counting gaps, missing data, or other ambiguities)
#'     
#'     11-14) Columns containing exactly three specific  bases (ACG, ACT, AGT, CGT, again, not counting gaps, missing data, or other ambiguities)
#'     
#'     15) Columns containing A, C, G, and T
#'     
#' The 7 site patterns are:
#' 
#'     1-4) Columns containing only one base (A, C, G, or T)
#'     
#'     5) Columns containing exactly two bases
#'     
#'     6) Columns containing exactly three bases
#'     
#'     7) Columns containing A, C, G, and T
#' 
#' @references 
#' Lewis, P. O. (2014). Posterior Predictive Bayesian Phylogenetic Model Selection. Systematic Biology, 63, 309-321.
#' 
#' @return The value t(y).


doLewisBins <- function(aln.mat,n.bins) {
  # recover()
  
  # For ease of grepping
  aln_mat <- tolower(aln.mat)
  
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
    warning("There are gap only columns in this alignment")
    
    aln_mat <- aln_mat[,-has.none]
    
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
  
  if (n.bins == 15) {
    ydot <- c(length(has.a),length(has.c),length(has.g),length(has.t),length(has.a.c),length(has.a.g),length(has.a.t),length(has.c.g),length(has.c.t),length(has.g.t),length(has.a.c.g),length(has.a.c.t),length(has.a.g.t),length(has.c.g.t),length(has.a.c.g.t))
  } else if (n.bins == 7) {
    has.three <- c(has.a.c.g,has.a.c.t,has.a.g.t,has.c.g.t)
    has.two <- c(has.a.c,has.a.g,has.a.t,has.c.g,has.c.t,has.g.t)
    ydot <- c(length(has.a),length(has.c),length(has.g),length(has.t),length(has.two),length(has.three),length(has.a.c.g.t))
  } else {
    stop("n.bins only defined for n=7 or n=15")
  }
  if ( sum(ydot) != NCOL ) {
    stop("sets do not sum to entire alignment, serious error has occured somewhere")
  }
  return(ydot)
}


#' Equation 11 from Lewis (2014).
#'
#' @param pattern.counts     Numeric. A vector of site patterns, such as one computed by doLewisBins.
#' 
#' @details This function calculates t(y) as in Lewis (2014) equation 11, where y is a vector of pattern counts.
#' 
#' @references 
#' Lewis, P. O. (2014). Posterior Predictive Bayesian Phylogenetic Model Selection. Systematic Biology, 63, 309-321.
#' 
#' @return The value t(y).

fn.t <- function(pattern.counts) {
  # recover()
  ns <- sum(pattern.counts) # number of sites in each alignment, nomenclature following Lewis 2014 for clarity
  return(-log(ns) + 1/ns * sum(pattern.counts*log(pattern.counts))) # equation 11 from Lewis 2014
}

#' Calculate the posterior predictive test of model fit described in Lewis (2014).
#'
#' @param aln.in     Character. An alignment in matrix form, as read in with readAlignment().
#' @param pp.aln     List. A list of (posterior predictive) alignments in matrix form (readAlignment() provides an alignment in matrix format).
#' @param phi     Numeric. Determines the relative importance of goodness-of-fit versus variance in the predicitve distribution, 1 specifies equal weight. Default 1.
#' @param bins     Numeric. How many site patterns should be considered, 7 or 15? Default 15.
#' 
#' @details This function calculates Lewis' GG (Gelfand-Ghosh) posterior test for model adequacy.
#' The test is designed to account for both the fit of the model to the data and the variance of the predicitive distribution of the model.
#' The argument phi weights goodness-of-fit against variance, with one specifying equal weights and increasing values putting more weight on goodness-of-fit.
#' The function is designed to take the real data and the simulated data and provide the value of the GG criterion.
#' To use this function to compare models, calculate the GG criterion for each model under consideration, choose the model with the minimum value.
#' 
#' Site patterns are binned so as to avoid issues when taking the log of the number of occurences (see Lewis equation 11).
#' The 15 site patterns are:
#' 
#'     1-4) Columns containing only one base (A, C, G, or T), not counting gaps, missing data, or other ambiguities
#'     
#'     5-10) Columns containing exactly two specific bases (AG, CT, etc, again, not counting gaps, missing data, or other ambiguities)
#'     
#'     11-14) Columns containing exactly three specific  bases (ACG, ACT, AGT, CGT, again, not counting gaps, missing data, or other ambiguities)
#'     
#'     15) Columns containing A, C, G, and T
#'     
#' The 7 site patterns are:
#' 
#'     1-4) Columns containing only one base (A, C, G, or T)
#'     
#'     5) Columns containing exactly two bases
#'     
#'     6) Columns containing exactly three bases
#'     
#'     7) Columns containing A, C, G, and T
#' 
#' @references 
#' Lewis, P. O. (2014). Posterior Predictive Bayesian Phylogenetic Model Selection. Systematic Biology, 63, 309-321.
#' 
#' @return The value of Lewis' GG criterion.

doLewisGelfandGhoshAdequacyTest <- function(real.aln,sim.alns,phi=1,bins=15) {
  # recover()
  nc <- bins # nc could be the number of site patterns possible, but since we consider only groupings, it is the number of groupings/bins
  ns <- ncol(real.aln) # number of sites in each alignment, nomenclature following Lewis 2014 for clarity
  np <- length(sim.alns) # number of simulated datasets, nomenclature following Lewis 2014 for clarity
  sim.counts <- matrix(ncol=nc,nrow=np) # a matrix with site pattern counts from simulated data
  y <- doLewisBins(real.aln,bins) # the vector of site pattern counts for the real data
  for (i in 1:np) {
    sim.counts[i,] <- doLewisBins(sim.alns[[i]],bins)
  }
  mu <- colSums(sim.counts)/np # mu is the vector of the average number of occurrences for the site patterns, nomenclature following Lewis 2014 for clarity
  GGp <- 2*ns*( (1/np * sum(apply(sim.counts,1,fn.t))) - fn.t(mu) ) # Equation 9
  GGg <- 2*ns*(1+phi)*( (fn.t(mu)+phi*fn.t(y))/(phi+1) - fn.t((mu+phi*y)/(phi+1)) ) # Equation 10
  return(GGp+GGg) # Equation 8
}

#' Calculate the posterior predictive test of model fit described in Lewis (2014).
#'
#' @param aln.in     Character. An alignment in matrix form, as read in with readAlignment().
#' @param sim.counts     Character matrix. The counts from all simulated alignments.
#' @param phi     Numeric. Determines the relative importance of goodness-of-fit versus variance in the predicitve distribution, 1 specifies equal weight. Default 1.
#' 
#' @details This function calculates the value of Lewis' Gelfand-Ghosh criterion given that the site patterns have been counted already with calculateNullDistributions.
#' Does not need number of bins as argument (this is defined by the number of site patterns already calculated)
#' 
#' @references 
#' Lewis, P. O. (2014). Posterior Predictive Bayesian Phylogenetic Model Selection. Systematic Biology, 63, 309-321.
#' 
#' @return The value of Lewis' GG criterion.

calculateGelfandGhoshCriterion <- function(real.counts,sim.counts,phi=1) {
  # recover()
  nc <- ncol(sim.counts) # nc could be the number of site patterns possible, but since we consider only groupings, it is the number of groupings/bins
  ns <- sum(sim.counts[1,]) # number of sites in each alignment, nomenclature following Lewis 2014 for clarity
  np <- nrow(sim.counts) # number of simulated datasets, nomenclature following Lewis 2014 for clarity
  y <- real.counts # the vector of site pattern counts for the real data
  if ( any(y == 0) ) {
    stop("Please reduce number of bins and try again. To go from 15 to 7, combine bins 5-10 and 11-14.")
  }
  mu <- colSums(sim.counts)/np # mu is the vector of the average number of occurrences for the site patterns, nomenclature following Lewis 2014 for clarity
  GGp <- 2*ns*( (1/np * sum(apply(sim.counts,1,fn.t))) - fn.t(mu) ) # Equation 9
  GGg <- 2*ns*(1+phi)*( (fn.t(mu)+phi*fn.t(y))/(phi+1) - fn.t((mu+phi*y)/(phi+1)) ) # Equation 10
  return(GGp+GGg) # Equation 8
}

# Internal function to calculate the number of times each nucleotide occurs at each site, averaging across uncertainties.

getCountsDistributeUncertainties <- function(vec) {
  GAP <- length(grep("-",vec,fixed=T)) + length(grep("?",vec,fixed=T))
  .X <- length(grep("x",vec,fixed=T))
  .N <- length(grep("n",vec,fixed=T))
  .M <- length(grep("m",vec,fixed=T))
  .R <- length(grep("r",vec,fixed=T))
  .W <- length(grep("w",vec,fixed=T))
  .S <- length(grep("s",vec,fixed=T))
  .Y <- length(grep("y",vec,fixed=T))
  .H <- length(grep("h",vec,fixed=T))
  .D <- length(grep("d",vec,fixed=T))
  .B <- length(grep("b",vec,fixed=T))
  .V <- length(grep("v",vec,fixed=T))
  .K <- length(grep("k",vec,fixed=T))
  
  nA <- length(grep("a",vec,fixed=T)) + 0.25 * GAP + 0.25 * .X + 0.25 * .N + 0.5 * .M + 0.5 * .R + 0.5 * .W + 0.33 * .V + 0.33 * .H + 0.33 * .D
  nC <- length(grep("c",vec,fixed=T)) + 0.25 * GAP + 0.25 * .X + 0.25 * .N + 0.5 * .M + 0.5 * .S + 0.5 * .Y + 0.33 * .V + 0.33 * .H + 0.33 * .B
  nG <- length(grep("g",vec,fixed=T)) + 0.25 * GAP + 0.25 * .X + 0.25 * .N + 0.5 * .R + 0.5 * .S + 0.5 * .K + 0.33 * .V + 0.33 * .D + 0.33 * .B
  nT <- length(grep("t",vec,fixed=T)) + 0.25 * GAP + 0.25 * .X + 0.25 * .N + 0.5 * .W + 0.5 * .Y + 0.5 * .K + 0.33 * .H + 0.33 * .D + 0.33 * .B
  return(c(nA,nC,nG,nT))
}

#' Computes values of given test statistics for each alignment in a set of posterior predictive alignments.
#' 
#' @param pps.alns     Character or List. Either the (absolute or relative) path to a directory containing posterior predictive alignments from Rev, or a list of alignments as character matrices.
#' 
#' @param statistics     Function. List of functions to be applied to each posterior predictive alignment.
#'
#' @param output.prefix     Character (default NULL). Either an overall prefix to be used for file names when writing distributions, or one prefix per statistic calculated. If NULL, no prefixes.
#' 
#' @param write     Boolean (default TRUE). Should the distributions of statistics be written to files? If FALSE, they are returned.
#' 
#' @param write.to     Character (default NA). Where to write statistics. Defaults to directory containing pps.alns (if reading alignments) or working directory (if not reading alignments).
#' 
#' @param store.alignments     Boolean (default FALSE). Should alignments be returned at end of function call?
#' 
#' @details Files are written to the directory that contains pps.dir.
#' If the user does not specify a unique set of output file prefixes, files are numbered in the order that the statistics are passed to calculateNullDistributions.
#' 
#' There are four currently implemented statistics that can be implemented
#'
#'  1) Bollback's unconstrained likelihood (unconstrainedLikelihood)
#'  
#'  2) Huelsenbeck et al.'s chi squared homogeneity (chiSquaredHomogeneity)
#'  
#'  3) Lartillot et al.'s biochemical diversity (biochemicalDiversity)
#'  
#'  4) 
#'  
#' Additionally, it is possible to use this method to simultaneously calculate the bin patterns necessary to evaluate using Lewis' Gelfand-Ghosh criterion.
#' In this case, make a wrapper function for doLewisBins and feed that function to calculateNullDistributions.
#' Then use calculateLewisGG to complete the test.
#' 
#' @examples
#' 
#' ## Make a function for doing Lewis' binning on the data
#' fnDoLewisBins7 <- function(aln.mat) {
#'  return(doLewisBins(aln.mat,7))
#' }
#' 
#' ## Devise some test statistic to calculate
#' fnNumSites <- function(aln.mat) {
#'  n.sites <- ncol(aln.mat)
#'  return(n.sites)
#' }
#' 
#' statistics <- list(unconstrainedLikelihood,fnDoLewisBins7,fnNumSites)
#' 
#' prefixes <- c("unconstrained_likelihood","lewis_bins_7","number_of_sites")
#'
#' ## Calculate the distributions and print to files, don't store alignments 
#' calculateNullDistribution(pps.alns="~/path/to/pps/alignments",statistics=statistics,output.prefix=prefixes)
#' 
#' @references 
#' Goldman, N. (1993). Statistical Tests of Models of DNA Substitution. Journal of Molecular Evolution, 36, 182-198.
#' 
#' Bollback, J. P. (2002). Bayesian Model Adequacy and Choice in Phylogenetics. Molecular Biology and Evolution, 19, 1171-1182.
#'
#' Huelsenbeck, J. P., Ronquist, F., Nielsen, R., and Bollback, J. P. (2001). Bayesian Inference of Phylogeny and Its Impact on Evolutionary Biology. Science, 294, 2310-2314.
#' 
#' Lartillot, N., Brinkmann, H., and Philippe, H. (2007). Suppression of long-branch attraction artefacts in the animal phylogeny using a site-heterogeneous model. BMC Evolutionary Biology, 7, 1-14.
#' 
#' Lewis, P. O. (2014). Posterior Predictive Bayesian Phylogenetic Model Selection. Systematic Biology, 63, 309-321.
#' 
#' @return If write = FALSE, the distributions of the test statistics. 
#' If store.alignments = TRUE, also (or solely, if write = TRUE) a list of the posterior predictive alignments.

calculateNullDistributions <- function(pps.alns,
                                       statistics,
                                       output.prefix,
                                       write = TRUE,
                                       write.to = NA,
                                       store.alignments = FALSE) {       
  
  ############################################
  # Set up for reading/processing alignments #
  ############################################
  # recover()
  # Check if we have alignments or folder in pps.aln
  if ( class(pps.alns) == "list" && class(pps.alns[[1]]) == "matrix" ) {
    read <- FALSE
    # Make a list (of lists) for storing the values we calculate on each alignment
    # One list per test statistic, one list element in that list per value calculated on simulated alignment
    simulated_values <- vector("list",length(statistics))
    for (i in 1:length(statistics)) {
      simulated_values[[i]] <- vector("list",length(pps.alns))
    }
    
  } else {
    if ( !file.exists(pps.alns) ) {
      stop(paste("There is no directory named ",pps.alns,sep=""))
    }
    # Get all folders with alignments
    read <- TRUE
    simulations <- list.files(pps.alns,full.names=TRUE)
    if (store.alignments) {
      # Make a list for alignments to go in
      alignments <- as.list(1:length(simulations))
    }
    # Make a list (of lists) for storing the values we calculate on each alignment
    # One list per test statistic, one list element in that list per value calculated on simulated alignment
    simulated_values <- vector("list",length(statistics))
    for (i in 1:length(statistics)) {
      simulated_values[[i]] <- vector("list",length(simulations))
    }
  }
  
  ################################################################
  # Gather values of test statistics from each simulated dataset #
  ################################################################
  # recover()
  # Set up progress bar
  pb <- txtProgressBar(min=0,max=length(simulations),style=3)
  if ( read ) {
    for (i in 1:length(simulations)) {
      
      # Find all parts of the alignment and read them in
      replicate_components <- list.files(simulations[i],full.names=TRUE)
      replicate_components_alns <- lapply(replicate_components,readAlignment)
      
      # Combine partitions into an overall alignment
      aln_mat <- do.call(cbind,replicate_components_alns)
      
      # Calculate test statistics
      for (j in 1:length(statistics)) {
        # simulated_values[[j]][[i]] <- do.call(statistics[j],list(aln_mat))
        simulated_values[[j]][[i]] <- statistics[[j]](aln_mat)
      }
      
      # Store the alignment if user requested
      if (store.alignments) {
        alignments[[i]] <- aln_mat
      }
      
      setTxtProgressBar(pb,i)
      
    }
  } else {
    for (i in 1:length(simulated_values[[1]])) {
      
      aln_mat <- pps.alns[[i]]
      
      # Calculate test statistics
      for (j in 1:length(statistics)) {
        # simulated_values[[j]][[i]] <- do.call(statistics[j],list(aln_mat))
        simulated_values[[j]][[i]] <- statistics[[j]](aln_mat)
      }
    }
    # Store the alignment if user requested
    if (store.alignments) {
      alignments <- pps.alns
    }
    
    setTxtProgressBar(pb,i)
    
  }
  
  close(pb)
  
  ###################################
  # Return or write requested items #
  ###################################
  
  # Ensure that we have a sufficient number of unique file name prefixes for writing
  if ( length(output.prefix) == 1 ) {
    if ( length(statistics) > 1 ) {
      output.prefix <- paste(output.prefix,1:length(statistics),sep="")
    }
  } else if ( ! length(unique(output.prefix)) == length(statistics) ) {
    warning("Not enough unique output prefixes specified, ignoring input, numbering files.")
    output.prefix <-1:length(statistics)
  }
  
  # Determine where to write distributions
  if ( !is.na(write.to) ) {
    out.files <- paste(write.to,"/",output.prefix,".txt",sep="")
  } else {
    if ( ! read ) {
      out.files <- paste(getwd(),"/",output.prefix,".txt",sep="")
    } else {
      out.files <- paste(dirname(pps.alns),"/",output.prefix,".txt",sep="")
    }
  }
  
  # recover()
  if ( store.alignments ) {
    if ( ! write ) {
      return(list(alignments=alignments,distributions=simulated_values))
    } else {
      for (i in 1:length(statistics)) {
        mat <- do.call(rbind,simulated_values[[i]])
        write.table(mat,file=out.files[i],sep=",",row.names=FALSE,col.names=FALSE,quote=FALSE)
      }
      return(alignments)
    }
  } else {
    if ( ! write ) {
      return(simulated_values)
    } else {
      for (i in 1:length(statistics)) {
        mat <- do.call(rbind,simulated_values[[i]])
        write.table(mat,file=out.files[i],sep=",",row.names=FALSE,col.names=FALSE,quote=FALSE)
      }
    }
  }
}
