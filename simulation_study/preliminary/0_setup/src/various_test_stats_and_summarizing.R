allMI <- function(phydat,exclude.invariant=TRUE) {
  # recover()
  
  # For ease of grepping, we make the alignment columns into lowercase strings
  aln_mat <- tolower(as.character(phydat))
  
  if (exclude.invariant) {
    aln_mat <- aln_mat[,!(lengths(apply(aln_mat,2,table)) == 1)]
  }
  
  A <- c("a","c","g","t")
  AA <- c("aa","ac","ag","at",
          "ca","cc","cg","ct",
          "ga","gc","gg","gt",
          "ta","tc","tg","tt")
  n <- dim(aln_mat)[1]
  
  # All sites' independent distributions
  p_x <- apply(aln_mat,2,function(x){
    sapply(A,function(a){
      sum(x == a)/n
    })
  })
  
  # joint distributions and MI
  mi <- sapply(1:(dim(aln_mat)[2]-1),function(i){
    sapply((i+1):dim(aln_mat)[2],function(j){
      # Exploit matrix math to get joint distribution
      # Makes a matrix where the rows are state in aln site j, cols in site i
      q_xy <- p_x[,j] %*% t(p_x[,i])
      
      # Flattens by column, so this puts aa, ac, ag, ...
      dim(q_xy) <- NULL
      
      # This produces p_xy in the same column vs row orientation as q
      p_xy <- sapply(A,function(ai){
        sapply(A,function(aj){
          sum((aln_mat[,i] == ai) & (aln_mat[,j] == aj))/n
        })
      })

      dim(p_xy) <- NULL
      
      # ij <- apply(aln_mat[,c(i,j)],1,paste0,collapse="")
      # p_xy <- sapply(AA,function(aa){
      #   sum(ij == aa)/n
      # })
      
      I_xy <- sum((p_xy * log(p_xy/q_xy))[!(p_xy == 0)])
      
      
      return(I_xy)
      
    })
  })
  return(unlist(mi))
}

TajimaD <- function(phydat) {
  # recover()
  
  # For ease of grepping, we make the alignment columns into lowercase strings
  aln_mat <- tolower(as.character(phydat))
  
  n <- dim(aln_mat)[1]
  
  n_inv <- sum(lengths(apply(aln_mat,2,table)) == 1)
  S <- dim(aln_mat)[2] - n_inv
  
  pair_d <- numeric(choose(n,2))
  index <- 0
  for (i in 1:(n - 1)) {
    for (j in (i+1):n) {
      index <- index + 1
      pair_d[index] <- sum(aln_mat[i,] != aln_mat[j,])
    }
  }
  
  PI <- sum(pair_d)/choose(n,2)
  
  a1 <- (sum(1/(1:(n - 1))))
  a2 <- (sum(1/((1:(n - 1)^2))))
  
  num <- PI - S/a1
  
  e1 <- 1/a1 * ((n+1)/(3*(n-1)) - 1/a1)
  e2 <- 1/(a1^2 + a2) * ( ((2*(n^2 + n + 3)) / (9*n*(n-1))) - ((n + 2)/(n*a1)) + (a2/(a1^2)) )
  
  denom <- sqrt(e1 * S + e2 * S * (S - 1))
  
  return(num/denom)
  
}

avgProportionSharedSegSites <- function(phydat) {
  aln <- tolower(as.character(phydat))
  n <- dim(aln)[1]
  
  tab <- apply(aln,2,table)
  
  S <- sum(lengths(tab) > 1)
  
  aln <- aln[,lengths(tab) > 1]
  
  prop_shared_pairwise <- sum(sapply(1:(n-1),function(i){
    sum(sapply((i+1):n,function(j){
      sum(aln[i,] == aln[j,])/S
    }))
  }))
  
  return(prop_shared_pairwise/choose(n,2))
  
}

goldman <- function(phydat) {
  N_xi <- attributes(phydat)$weight
  N <- sum(N_xi)
  sum(N_xi * log(N_xi)) - N * log(N)
}

SEV <- function(phydat) {
  aln <- tolower(as.character(phydat))
  n <- dim(aln)[1]
  ent <- apply(aln,2,function(x){
    p <- table(x)/n
    return(sum(p * log(p)))
  })
  return(var(ent))
}

SES <- function(phydat) {
  aln <- tolower(as.character(phydat))
  n <- dim(aln)[1]
  ent <- apply(aln,2,function(x){
    p <- table(x)/n
    return(sum(p * log(p)))
  })
  return(skew(ent))
}

allSiteEntropies <- function(phydat) {
  aln <- tolower(as.character(phydat))
  n <- dim(aln)[1]
  ent <- apply(aln,2,function(x){
    p <- table(x)/n
    return(sum(p * log(p)))
  })
  return(ent)
}

eSFS <- function(phydat) {
  aln <- tolower(as.character(phydat))
  n <- dim(aln)[1]
  min_site_count <- apply(aln,2,function(x){
    return(min(table(x)))
  })
  esfs <- sapply(1:n,function(i){sum(min_site_count == i)})
  return(esfs)
}

bowker <- function(seq1,seq2,return.p.val=TRUE) {
  nt <- c("a","c","g","t")
  
  # # We don't need this matrix, we can calculate N[i,j] at the same time we calculate the test statistic
  # N <- matrix(nrow=4,ncol=4)
  # for (i in 1:4) {
  #   for (j in 1:4) {
  #     N[i,j] <- sum(seq1 == nt[i] & seq2 == nt[j])
  #   }
  # }
  
  test_stat <- sum(sapply(1:3,function(i){
    sum(sapply((i+1):4,function(j){
      n_ij <- seq1 == nt[i] & seq2 == nt[j]
      n_ji <- seq1 == nt[j] & seq2 == nt[i]
      return( ((n_ij - n_ji)^2)/(n_ij + n_ji) )
    }))
  }))
  
  if (return.p.val) {
    return(pchisq(test_stat,df=????))
  } else {
    return(test_stat)
  }
}

stuart <- function(seq1,seq2,return.p.val=TRUE) {
  nt <- c("a","c","g","t")
  
  # We don't need this matrix, we can calculate N[i,j] at the same time we calculate the test statistic
  N <- matrix(nrow=4,ncol=4)
  for (i in 1:4) {
    for (j in 1:4) {
      N[i,j] <- sum(seq1 == nt[i] & seq2 == nt[j])
    }
  }
  
  V <- matrix(nrow=3,ncol=3)
  for (i in 1:3) {
    V[i,j] <- 
    for (j in c(1:3)[!i]) {
      V[i,j] <- 
    }
  }
  
  if (return.p.val) {
    return(pchisq(test_stat,df=????))
  } else {
    return(test_stat)
  }
}

SRHChiSquaredMarginal <- function(phydat) {
  aln <- tolower(as.character(phydat))
  n <- dim(aln)[1]
  
  all_comps <- numeric(length=choose(n,2))
  index <- 0
  
  # Calculate all stats
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      index <- index + 1
      all_comps[index] <- bowker(aln[i,],aln[j,])
    }
  }
  
  return(???)
  
}

kurt <- function(x) {
  mu <- mean(x)
  n <- length(x)
  ku <- (sum((x-mu)^4)/n) / (sum((x-mu)^2)/n)^2
  return(ku)
}

# Internal function to calculate the skewness of a distribution.

skew <- function(x) {
  mu <- mean(x)
  n <- length(x)
  sk <- (sum((x-mu)^3)/n) / (sum((x-mu)^2)/n)^3/2
  return(sk)
}

dir1 <- "simulation_study/preliminary/1_decreasing_site_counts/frac-0/data"
dir2 <- "simulation_study/preliminary/2_epistasis/prop-0.5_d-0/data"
dir3 <- "simulation_study/preliminary/2_epistasis/prop-0.5_d-10/data"

true.tree <- read.tree("data/flu_HA.tre")

aln1 <- lapply(list.files(dir1,full.names=TRUE),function(aln){
  read.phyDat(aln,"nexus")
})

aln2 <- lapply(list.files(dir2,full.names=TRUE),function(aln){
  read.phyDat(aln,"nexus")
})

aln3 <- lapply(list.files(dir3,full.names=TRUE),function(aln){
  read.phyDat(aln,"nexus")
})


nulldist <- unlist(lapply(aln1,function(aln){
  fitch(true.tree,aln)
}))

epidist <- unlist(lapply(aln2,function(aln){
  fitch(true.tree,aln)
}))

epidist2 <- unlist(lapply(aln3,function(aln){
  fitch(true.tree,aln)
}))


range(nulldist)
range(epidist)


allbeans <- list(non_epistatic=nulldist,
                 prop_0.5_d_0=epidist,
                 prop.0.5_d_10=epidist2)
beanplot(allbeans,cutmin=0,what=c(0,1,1,0),col=list("#66666690","#FF000090","#0000FF90"),ylab="Fitch score of alignment")



reg.lb <- lapply(aln1,function(aln){
  doLewisBins(as.character(aln),15)
})
reg.lb <- do.call(rbind,reg.lb)

epi.lb <- lapply(aln2,function(aln){
  doLewisBins(as.character(aln),15)
})
epi.lb <- do.call(rbind,epi.lb)

epi2.lb <- lapply(aln3,function(aln){
  doLewisBins(as.character(aln),15)
})
epi2.lb <- do.call(rbind,epi2.lb)


colMeans(reg.lb)
colMeans(epi.lb)
colMeans(epi2.lb)


rlb <- vector("list",15)
for (i in 1:15) {
  rlb[[i]] <- reg.lb[,i]
}

elb <- vector("list",15)
for (i in 1:15) {
  elb[[i]] <- epi.lb[,i]
}

elb2 <- vector("list",15)
for (i in 1:15) {
  elb2[[i]] <- epi2.lb[,i]
}


reg.reduced.lb <- cbind(rowSums(reg.lb[,1:4]),rowSums(reg.lb[,5:10]),rowSums(reg.lb[,11:14]),reg.lb[,15])
epi.reduced.lb <- cbind(rowSums(epi.lb[,1:4]),rowSums(epi.lb[,5:10]),rowSums(epi.lb[,11:14]),epi.lb[,15])
epi2.reduced.lb <- cbind(rowSums(epi2.lb[,1:4]),rowSums(epi2.lb[,5:10]),rowSums(epi2.lb[,11:14]),epi2.lb[,15])

colMeans(reg.reduced.lb)
colMeans(epi.reduced.lb)
colMeans(epi2.reduced.lb)

rrlb <- vector("list",4)
for (i in 1:4) {
  rrlb[[i]] <- reg.reduced.lb[,i]
}

erlb <- vector("list",4)
for (i in 1:4) {
  erlb[[i]] <- epi.reduced.lb[,i]
}

e2rlb <- vector("list",4)
for (i in 1:4) {
  e2rlb[[i]] <- epi2.reduced.lb[,i]
}


allrbeans <- vector("list",12)
for (i in 1:4) {
  allrbeans[[3*i-2]] <- rrlb[[i]]
  allrbeans[[3*i-1]]   <- erlb[[i]]
  allrbeans[[3*i]]   <- e2rlb[[i]]
}
names(allrbeans) <- c(1,1,1,2,2,2,3,3,3,4,4,4)

beanplot(allrbeans,cutmin=0,what=c(0,1,1,0),col=list("#66666690","#FF000090","#0000FF90"),xlab="# distinct nucleotides at a site")
legend("topright",legend=c("all IID sites","50% epistatic, d=0","50% epistatic, d=10"),fill=c("#66666690","#FF000090","#0000FF90"),bty="n",border=NA)




reg.sfs <- lapply(aln2,function(aln){
  eSFS(aln)
})
reg.sfs <- do.call(rbind,reg.sfs)

epi.sfs <- lapply(aln1,function(aln){
  eSFS(aln)
})
epi.sfs <- do.call(rbind,epi.sfs)


colMeans(reg.sfs)
colMeans(epi.sfs)

rsfs <- vector("list",25)
for (i in 1:24) {
  rsfs[[i]] <- reg.sfs[,i]
}
rsfs[[25]] <- reg.sfs[,49]

esfs <- vector("list",25)
for (i in 1:24) {
  esfs[[i]] <- epi.sfs[,i]
}
esfs[[25]] <- epi.sfs[,49]

allsfs <- vector("list",50)
for (i in 1:25) {
  allsfs[[2*i-1]] <- rsfs[[i]]
  allsfs[[2*i]]   <- esfs[[i]]
}

beanplot(allsfs[-c(1:2,49:50)],cutmin=0,what=c(0,1,1,0),col=list("#66666690","#FF000090"))

hist(reg.sfs[,2],col="#66666690")
hist(epi.sfs[,2],col="#FF000090",add=TRUE)

pdf("~/Pictures/sfs_comparison.pdf",width=8,height=4)
  par(mfrow=c(1,2))
  beanplot(rsfs,cutmin=0,what=c(0,1,1,0),col=list("#66666690"))
  beanplot(esfs,cutmin=0,what=c(0,1,1,0),col=list("#FF000090"))
dev.off()

