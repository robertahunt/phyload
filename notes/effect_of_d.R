library(expm)
library(phangorn)

########
# Functions we need
########

# Make the instantaneous rate matrix for the epistatic model on site-pairs
assembleEpiQ <- function(er,df,epistasis_d) {
  
  S <- matrix(0,4,4)
  S[lower.tri(S)] <- er
  S <- t(S)
  S[lower.tri(S)] <- er
  
  d1 <- c(1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4)
  d2 <- c(1,2,3,4,1,2,3,4,1,2,3,4,1,2,3,4)
  
  unscaled_Q <- matrix(0,16,16)
  mu <- 0
  for (i in 1:16) {
    # x_1 and x_2 help us tell what cell in the GTR model we'd be in
    # Namely, they tell us the first nucleotide and second nucleotide in the "from" doublet under consideration
    x_1 = d1[i]
    x_2 = d2[i]
    
    for (j in 1:16) {
      # y_1 and y_2 tell us the same thing for the "to" doublet
      y_1 = d1[j]
      y_2 = d2[j]
      
      if (i == j) { # Catch diagonal entries first, this allows us to not add exceptions to our cases in next if statements
        unscaled_Q[i,j] = 0
      } else if ( (x_1 == 1 && x_2 == 4 || x_1 == 4 && x_2 == 1 || x_1 == 2 && x_2 == 3 || x_1 == 3 && x_2 == 2) && (y_1 == 1 && y_2 == 4 || y_1 == 4 && y_2 == 1 || y_1 == 2 && y_2 == 3 || y_1 == 3 && y_2 == 2) ) {  # Change from one canonically paired doublet to another
        unscaled_Q[i,j] <- abs(epistasis_d * S[x_1,y_1] * S[x_2,y_2] * df[j])
        mu <- mu + df[i] * unscaled_Q[i,j]
      } else if (x_2 == y_2) { # single base change at first base, second base is the same
        unscaled_Q[i,j] <- abs(S[x_1,y_1] * df[j])
        mu <- mu + df[i] * unscaled_Q[i,j] * 0.5
      } else if (x_1 == y_1) { # single base change at second base, first base is the same
        unscaled_Q[i,j] <- abs(S[x_2,y_2] * df[j])
        mu <- mu + df[i] * unscaled_Q[i,j] * 0.5
      } else { # double mutation of a disallowed variety
        unscaled_Q[i,j] = 0
      }
    }
  }
  
  Q_epi <- 1/mu * unscaled_Q
  
  diag(Q_epi) <- -rowSums(Q_epi)
  
  return(Q_epi)
}

# Calculate the proportion of substitutions that are doublet substitutions
rateProportionDoublets <- function(er,df,epistasis_d) {
  # recover()
  
  Q_epi <- assembleEpiQ(er,df,epistasis_d)
  
  sum_doublet <- 0
  sum_single <- 0
  
  d1 <- c(1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4)
  d2 <- c(1,2,3,4,1,2,3,4,1,2,3,4,1,2,3,4)
  
  for (i in 1:16) {
    # x_1 and x_2 help us tell what cell in the GTR model we'd be in
    # Namely, they tell us the first nucleotide and second nucleotide in the "from" doublet under consideration
    x_1 = d1[i]
    x_2 = d2[i]
    
    for (j in 1:16) {
      # y_1 and y_2 tell us the same thing for the "to" doublet
      y_1 = d1[j]
      y_2 = d2[j]
      if (i != j) {
        if ( (x_1 == 1 && x_2 == 4 || x_1 == 4 && x_2 == 1 || x_1 == 2 && x_2 == 3 || x_1 == 3 && x_2 == 2) && (y_1 == 1 && y_2 == 4 || y_1 == 4 && y_2 == 1 || y_1 == 2 && y_2 == 3 || y_1 == 3 && y_2 == 2) ) {  # Change from one canonically paired doublet to another
          sum_doublet <- sum_doublet + Q_epi[i,j]
        } else if (x_2 == y_2) { # single base change at first base, second base is the same
          sum_single <- sum_single + Q_epi[i,j]
        } else if (x_1 == y_1) { # single base change at second base, first base is the same
          sum_single <- sum_single + Q_epi[i,j]
        }
      }
    }
  }
  return(sum_doublet/(sum_doublet + sum_single))
}

# # Make a 13x13 matrix out of the 16x16 matrix, where the new state is A at doublet site 1, or G at doublet site 2
# # Will return matrix with condensed state in first row/column
# condenseEpiTransitionProbabilityMatrix <- function(Q_epi,acgt,first_or_second) {
#   # recover()
#   
#   d1 <- c(1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4)
#   d2 <- c(1,2,3,4,1,2,3,4,1,2,3,4,1,2,3,4)
#   
#   if ( is.numeric(acgt) ) {
#     if (acgt < 1 || acgt > 4) {
#       stop("Invalid state to condense on.")
#     } else {
#       state <- acgt
#     }
#   } else {
#     acgt <- toupper(acgt)
#     if (acgt == "A") {
#       state <- 1
#     } else if (acgt == "C") {
#       state <- 2
#     } else if (acgt == "G") {
#       state <- 3
#     } else if (acgt == "T") {
#       state <- 4
#     } else {
#       stop("Invalid state to condense on.")
#     }
#   }
#   
#   if ( first_or_second == 1 ) {
#     condenser <- d1
#   } else if ( first_or_second == 2 ) {
#     condenser <- d2
#   } else {
#     stop("Option first_or_second takes a numeric argument")
#   }
#   
#   new_prob_to_state <- numeric(16)
#   for (i in 1:16) {
#     if ( condenser[i] == state ) {
#       new_prob_to_state[i] <- NA
#     } else {
#       new_prob_to_state[i] <- sum(Q_epi[i,condenser == state])
#     }
#   }
#   new_prob_to_state <- new_prob_to_state[!is.na(new_prob_to_state)]
#   
#   new_prob_from_state <- numeric(16)
#   for (j in 1:16) {
#     if ( condenser[j] == state ) {
#       new_prob_from_state[j] <- NA
#     } else {
#       new_prob_from_state[j] <- sum(Q_epi[condenser == state,j])
#     }
#   }
#   new_prob_from_state <- new_prob_from_state[!is.na(new_prob_from_state)]
#   
#   # Make new matrix
#   
#   # First, we start by dropping out all rows and columns that are going into the aggregate state
#   reduced_Q <- Q_epi[condenser != state,condenser != state]
#   
#   # Add in new rate of leaving
#   reduced_Q <- rbind(new_prob_from_state,reduced_Q)
#   
#   # Add in new rate of arriving
#   reduced_Q <- cbind(c(0,new_prob_to_state),reduced_Q)
#   
#   # Fix diagonal
#   diag(reduced_Q) <- -rowSums(reduced_Q)
#   
#   return(reduced_Q)
#   
# }

# expectedEpiPropInv <- function(er,df,epistasis_d,gamma_shape_parameter,branch_length) {
#   recover()
# 
#   d1 <- c(1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4)
#   d2 <- c(1,2,3,4,1,2,3,4,1,2,3,4,1,2,3,4)
#   
#   Q_epi <- assembleEpiQ(er,df,epistasis_d)
#   
#   if ( !is.numeric(gamma_shape_parameter)) {
#     P_epi <- expm::expm(Q_epi * branch_length)
#     
#     # To get nucleotide-level results, we must marginalize over each pair of the doublet in turn
#     # First, we consider the first site at a doublet, first site is A,C,G,T, and in all cases we marginalize over the second state
#     p_inv <- numeric(4)
#     
#     # First site
#     p_inv_1 <- numeric(4)
#     for (n in 1:4) {
#       total_stationary_frequency <- 0.5 * sum(df[d1 == n])
#       p_inv_1[n] <- total_stationary_frequency * sum(P_epi[d1 == n,d1 == n])
#     }
#     
#     # Second site
#     p_inv_2 <- numeric(4)
#     for (n in 1:4) {
#       total_stationary_frequency <- 0.5 * sum(df[d2 == n])
#       p_inv_1[n] <- total_stationary_frequency * sum(P_epi[d2 == n,d2 == n])
#     }
#     
#     p_inv <- (p_inv_1 + p_inv_2)/2
#     
#    return(p_inv)
#   } else {
# 
#     gamma_cats <- discrete.gamma(gamma_shape_parameter,4)
#     
#     p_inv <- matrix(0,4,4)
#     for (g in 1:4) {
#       P_epi <- expm::expm(Q_epi * branch_length * gamma_cats[g])
#       
#       # First site
#       p_inv_1 <- numeric(4)
#       for (n in 1:4) {
#         total_stationary_frequency <- 0.5 * sum(df[d1 == n])
#         p_inv_1[n] <- total_stationary_frequency * sum(P_epi[d1 == n,d1 == n])
#       }
#       
#       # Second site
#       p_inv_2 <- numeric(4)
#       for (n in 1:4) {
#         total_stationary_frequency <- 0.5 * sum(df[d2 == n])
#         p_inv_1[n] <- total_stationary_frequency * sum(P_epi[d2 == n,d2 == n])
#       }
#       
#       p_inv[g,] <- (p_inv_1 + p_inv_2)/2
#       
#     }
#     
#     return(colMeans(p_inv))
#       
#   }
#   
# }

# expectedEpiPropInv <- function(er,df,epistasis_d,gamma_shape_parameter,branch_length) {
#   recover()
#   
#   d1 <- c(1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4)
#   d2 <- c(1,2,3,4,1,2,3,4,1,2,3,4,1,2,3,4)
#   
#   Q_epi <- assembleEpiQ(er,df,epistasis_d)
#   
#   if ( !is.numeric(gamma_shape_parameter)) {
#     gamma_cats <- 1
#   } else {
#     gamma_cats <- discrete.gamma(gamma_shape_parameter,4)
#   }
#   
#   p_inv <- matrix(0,length(gamma_cats),4)
#   
#   for (g in 1:length(gamma_cats)) {
#     # First site
#     p_inv_1 <- numeric(4)
#     for (n in 1:4) {
#       submatrix <- Q_epi[d1 == n,d1 == n]
#       total_stationary_frequency <- 0.5 * sum(df[d1 == n])
#       total_leaving_rate <- -sum(submatrix)
#       p_inv_1[n] <- total_stationary_frequency * exp(-total_leaving_rate * branch_length * gamma_cats[g])
#     }
#     
#     # Second site
#     p_inv_2 <- numeric(4)
#     for (n in 1:4) {
#       submatrix <- Q_epi[d2 == n,d2 == n]
#       total_stationary_frequency <- 0.5 * sum(df[d2 == n])
#       total_leaving_rate <- sum(submatrix)
#       p_inv_2[n] <- total_stationary_frequency * exp(total_leaving_rate * branch_length * gamma_cats[g])
#     }
#     
#     p_inv[g,] <- (p_inv_1 + p_inv_2)
#     
#   }
#     
#   return(colMeans(p_inv))
# 
# }

expectedEpiPropInvDoublet <- function(er,df,epistasis_d,gamma_shape_parameter,branch_length) {
  # recover()
  
  d1 <- c(1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4)
  d2 <- c(1,2,3,4,1,2,3,4,1,2,3,4,1,2,3,4)
  
  Q_epi <- assembleEpiQ(er,df,epistasis_d)
  
  if ( !is.numeric(gamma_shape_parameter)) {
    gamma_cats <- 1
  } else {
    gamma_cats <- discrete.gamma(gamma_shape_parameter,4)
  }
  
  p_inv <- matrix(0,length(gamma_cats),16)
  
  for (g in 1:length(gamma_cats)) {
    P_epi <- expm::expm(Q_epi * branch_length * gamma_cats[g])
    p_inv[g,] <- df * diag(P_epi)
  }
  
  return(colMeans(p_inv))
  
}

expectedIIDPropInv <- function(er,bf,gamma_shape_parameter,branch_length) {
  
  # recover()
  
  # Assemble rate matrix
  Q <- matrix(0,4,4)
  Q[upper.tri(Q)] <- er
  Q <- t(Q)
  Q[upper.tri(Q)] <- er
  Q <- t(Q * bf)
  diag(Q) <- -rowSums(Q)
  mu <- -sum(bf * diag(Q))
  Q <- Q * 1/mu
  
  if ( class(gamma_shape_parameter) != "numeric") {
    gamma_cats <- 1
  } else {
    gamma_cats <- discrete.gamma(gamma_shape_parameter,4)
    
    p_inv <- matrix(NA,length(gamma_cats),4)
    
    for (g in 1:length(gamma_cats)) {
      # To be invariant, barring back-substitutions, we must start in a state and then never leave
      p_inv[g,] <- bf * exp(diag(Q*branch_length*gamma_cats[g]))
    }
    
    return(colMeans(p_inv))
  }
  
}

marginalizeDoubletFrequencies <- function(df) {
  d1 <- c(1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4)
  d2 <- c(1,2,3,4,1,2,3,4,1,2,3,4,1,2,3,4)
  bf <- numeric(4)
  for (i in 1:16) {
    bf[d1[i]] <- bf[d1[i]] + 0.5 * df[i] 
    bf[d2[i]] <- bf[d2[i]] + 0.5 * df[i] 
  }
  return(bf)
}

#######
# Calculate stuff!
#######

# Our flu exchange rates
er <- c(1.882161, 7.009179, 0.914813, 0.495852, 7.666181, 1.000000)
er <- er/sum(er)

# Flu alpha parameter for ASRV
gamma_alpha <- 0.440894

# Nasrallah and Huelsenbeck's exchange rates
S <- c(1,2,1,1,2,1)
S <- S/sum(S)
  
# Nasrallah and Huelsenbeck's doublet stationary frequencies
df <- c(0.015, 0.025, 0.011, 0.15, 0.0175, 0.019, 0.175, 0.014,0.014, 0.22, 0.01, 0.06, 0.14, 0.04, 0.065, 0.0245)
df <- df/sum(df)

# Tunicate parameters
er <- c(0.11,0.187,0.107,0.116,0.366,0.114)
er <- er/sum(er)

df <- c(0.01745,0.02179,0.0231,0.154,0.0192,0.02119,0.184,0.01484,0.01476,0.245,0.02152,0.05001,0.122,0.01751,0.05221,0.02157)
df <- df/sum(df)

bf <- c(0.319,0.209,0.245,0.227)
bf <- bf/sum(bf)

# Our values in the study
rateProportionDoublets(er,df,0.5)
rateProportionDoublets(er,df,2)
rateProportionDoublets(er,df,8)

# Large sequence
d <- exp(seq(log(1/1024),log(1024),length.out=100))

p <- sapply(d,function(this_d){rateProportionDoublets(er,df,this_d)})

plot(d,p,log="x")

# Invariant sites
expectedEpiPropInvDoublet(er,df,gamma_alpha,epistasis_d=0,branch_length=6)
expectedEpiPropInvDoublet(er,df,gamma_alpha,epistasis_d=0.5,branch_length=6)
expectedEpiPropInvDoublet(er,df,gamma_alpha,epistasis_d=2,branch_length=6)
expectedEpiPropInvDoublet(er,df,gamma_alpha,epistasis_d=8,branch_length=6)

epi_inv <- sapply(seq(0,100,1),function(d){
  expectedEpiPropInvDoublet(er,df,gamma_alpha,epistasis_d=d,branch_length=5.928407)
})

colSums(epi_inv)

plot(NULL,NULL,xlim=c(0,100),ylim=c(0,1),xlab="d",ylab="prop inv")

# lines(seq(0,100,1),epi_inv[1,],col="red")
# lines(seq(0,100,1),epi_inv[2,],col="green")
# lines(seq(0,100,1),epi_inv[3,],col="gold2")
# lines(seq(0,100,1),epi_inv[4,],col="blue")
# lines(seq(0,100,1),colSums(epi_inv),col="black")
# legend("topleft",fill=c("red","green","gold2","blue","black"),legend=c("A","C","G","T","total"),bty="n",border=NA)

plot(seq(0,100,1),colSums(epi_inv),type="l",ylim=c(0,1),xlab="d",ylab="prop inv")
abline(h=sum(expectedIIDPropInv(er,bf,gamma_alpha,5.928407)),col="grey")

