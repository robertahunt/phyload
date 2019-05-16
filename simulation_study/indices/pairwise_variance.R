#' Compute the posterior predictive test statistic of Goldman and Yang (1993)
#'
#' @param aln              phyDat alignment.
#' 
#' @details This function calculates the variance of the number of substitutions in all pairwise comparisons.
#' 
#' @return The value of the test statistic.

pairwiseVariance <- function(aln) {
  pairs <- lapply(1:(length(aln)-1),function(i){
    lapply((i+1):length(aln),function(j){
      sum(aln[[i]] != aln[[j]])
    })
  })
  return(var(unlist(pairs)))
}
