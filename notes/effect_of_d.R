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

# Tunicate parameters
er <- c(0.11,0.187,0.107,0.116,0.366,0.114)
er <- er/sum(er)

df <- c(0.01745,0.02179,0.0231,0.154,0.0192,0.02119,0.184,0.01484,0.01476,0.245,0.02152,0.05001,0.122,0.01751,0.05221,0.02157)
df <- df/sum(df)

bf <- c(0.319,0.209,0.245,0.227)
bf <- bf/sum(bf)

# Our values in the study
rateProportionDoublets(er,df,0.0)
rateProportionDoublets(er,df,0.5)
rateProportionDoublets(er,df,2)
rateProportionDoublets(er,df,8)
rateProportionDoublets(er,df,1000)

# Large sequence
d <- exp(seq(log(1/1000),log(10000),length.out=1000))

p <- sapply(d,function(this_d){rateProportionDoublets(er,df,this_d)})

# pdf("notes/figures/p_vs_d.pdf",)
#   par(lend=2)
  plot(d,p,type="l",log="x",ylab="proportion doublet")
# dev.off()
