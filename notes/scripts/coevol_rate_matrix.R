assembleCoEvolBF <- function(r,d,s,phi) {
  # Stationary frequencies appear to be simplex(nu) for
  # nu = ifelse(i in phi,d/s,1)
  bf <- rep(1,16)
  bf[bf %in% phi] <- d/s
  bf <- bf/sum(bf)
  return(bf)
}

assembleCoEvolMatrix <- function(r,d,s,phi) {
  Q_ce <- matrix(0,16,16)
  
  d1 <- c(1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4)
  d2 <- c(1,2,3,4,1,2,3,4,1,2,3,4,1,2,3,4)
  
  unscaled_Q <- matrix(0,16,16)
  mu <- 0
  for (i in 1:16) {
    for (j in 1:16) {
      if ( (i %in% phi) && !(j %in% phi) ) {
        # Leave pair
        Q_ce[i,j] <- s
      } else if ( !(i %in% phi) && (j %in% phi) ) {
        # Enter pair
        Q_ce[i,j] <- d
      } else if ( i != j ) {
        # Stay out of pair
        Q_ce[i,j] <- r
      }
    }
  }
  diag(Q_ce) <- -rowSums(Q_ce)
  
  bf <- assembleCoEvolBF(r,d,s,phi)
  
  mu <- -sum(diag(Q_ce) * bf)
  
  Q_ce <- Q_ce/mu
  
  return(Q_ce)
}

simCoEvol2SubsStartInPhi <- function(nsim,r,d,s,phi) {

  sapply(1:nsim,function(i){
    state_1 <- sample(phi,1)
    
    Q_ce <- assembleCoEvolMatrix(r,d,s,phi)
    
    rel_probs <- Q_ce[state_1,]
    rel_probs[state_1] <- 0.0
    
    state_2 <- sample(1:16,1,prob=rel_probs)
    
    rel_probs <- Q_ce[state_2,]
    rel_probs[state_2] <- 0.0
    
    state_3 <- sample(1:16,1,prob=rel_probs)
    
    return(state_3)
    
  })
}

sum(simCoEvol2SubsStartInPhi(100,0.5,100,1,c(4,7,10,13)) %in% c(4,7,10,13))


coev <- assembleCoEvolMatrix(0.5,100,1,c(4,7,10,13))

cbf <- assembleCoEvolBF(0.5,100,1,c(4,7,10,13))

P <- expm::expm(coev*0.1)


fwd_equals_bwd <- unlist(lapply(1:15,function(i){
  lapply((i+1):16,function(j){
    round(P[i,j] * cbf[i],6) == round(P[j,i] * cbf[j],6)
  })
}))



