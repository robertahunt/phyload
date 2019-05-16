#' Compute the posterior predictive test statistic of Goldman and Yang (1993)
#'
#' @param aln              phyDat alignment.
#' @param drop.ambiguous   TRUE/FALSE should columns with ambiguous characters be dropped?
#' 
#' @details This function calculates the value of Goldman's (1993) "unconstrained likelihood" to use as a test statistic, as in Bollback (2002).
#' Specifically, the value computed is given by equation 6b in Goldman (1993) and equation 7 in Bollback (2002).
#' 
#' @return The value of the test statistic.

unconstrainedLikelihood <- function(aln,drop.ambiguous=FALSE) {
  # Nomenclature follows Bollback's (2002) description of the test statistic.
  
  # For ease of grepping
  aln_mat <- tolower(as.character(aln))
  
  # Collapse sites into character strings, use table to compute number of occurences
  N_Xi <- table(apply(aln_mat,2,paste,collapse=""))
  
  if ( drop.ambiguous ) {
    has_ambiguity <- ceiling(which(!(N_Xi %in% c("a","c","g","t")))/length(N_Xi))
    N_Xi <- N_Xi[!has_ambiguity]
  }
  
  # N is the total number of sites
  N <- sum(N_Xi)
  
  # Compute statistic
  T_of_X <- sum(N_Xi * log(N_Xi)) - N*log(N)
  return(T_of_X)
}
