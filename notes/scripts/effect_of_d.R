library(expm)
library(phangorn)
library(viridis)

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

  
#######
# Compare how different the Q matrices are
#######
library(phylomd)

countPriorSubs <- function(Q,tree) {
  # recover()
  
  if ( class(Q) != "substitution.model" ) {
    stop("Q must be of class substitution.model")
  }
  
  nstates <- length(Q$states)
  
  counts <- numeric(choose(nstates,2))
  
  edges <- list(1:length(tree$edge.length))
  
  tips <- rep("-",length(tree$tip.label))
  
  idx <- 0
  pb <- txtProgressBar(style=3)
  for (i in 1:(nstates-1)) {
    for (j in (i+1):nstates) {
      idx <- idx + 1
      L <- matrix(0,nstates,nstates)
      L[i,j] <- 1
      L[j,i] <- 1
      counts[idx] <- phylo.nsubs.moments(tree,Q,L,edges,1,tips)[2]
      names(counts)[idx] <- paste0(Q$states[i],"<->",Q$states[j])
      setTxtProgressBar(pb,idx/length(counts))
    }
  }
  return(counts)
}

doubletCounts2siteCounts <- function(dc) {
  # recover()
  
  from <- do.call(rbind,strsplit(names(dc),"<->"))[,1]
  to <- do.call(rbind,strsplit(names(dc),"<->"))[,2]
  
  doublet_states <- unique(unlist(strsplit(names(dc),"<->")))
  single_states <- unique(unlist(strsplit(doublet_states,"")))
  
  counts <- matrix(0,length(single_states),length(single_states),dimnames=list(single_states,single_states))
  
  # get counts on single sites
  for (i in 1:length(dc)) {
    x <- strsplit(from[i],"")[[1]]
    y <- strsplit(to[i],"")[[1]]
    
    if ( x[1] == y[1] && x[2] != y[2] ) {
      counts[x[2],y[2]] <- counts[x[2],y[2]] + dc[i]
    } else if ( x[1] != y[1] && x[2] == y[2] ) {
      counts[x[1],y[1]] <- counts[x[1],y[1]] + dc[i]
    } else if ( x[1] != y[1] && x[2] != y[2] ) {
      counts[x[2],y[2]] <- counts[x[2],y[2]] + dc[i]
      counts[x[1],y[1]] <- counts[x[1],y[1]] + dc[i]
    }
  }
  
  # flatten counts
  flat_counts <- numeric(choose(length(single_states),2))
  idx <- 0
  for (i in 1:(length(single_states)-1)) {
    for (j in (i+1):length(single_states)) {
      idx <- idx + 1
      flat_counts[idx] <- counts[i,j] + counts[j,i]
      names(flat_counts)[idx] <- paste0(single_states[i],"<->",single_states[j])
    }
  }
  return(flat_counts)
}

# tunicate tree
tree <- read.tree("simulation_study/simulation_scripts/simulation_tree_tunicates.tre")

# root to 0-length branch, this is arbitrary but required
tree <- root(tree,outgroup=tree$tip.label[1],resolve.root=TRUE)

# What we would get if we simply expanded the GTR rate matrix into a 16x16 matrix
independent.df <- c(bf[1] * bf, bf[2] * bf, bf[3] * bf, bf[4] * bf)
bigger.gtr <- assembleEpiQ(er,independent.df,0)

d0.0 <- assembleEpiQ(er,df,0)
d0.5 <- assembleEpiQ(er,df,0.5)
d2.0 <- assembleEpiQ(er,df,2)
d8.0 <- assembleEpiQ(er,df,8)
d1000.0 <- assembleEpiQ(er,df,1000)

# make "substitution model" class objects for use in mapping
bigger.jc <- list(Q=assembleEpiQ(rep(1/6,6),rep(1/16,16),0),
                  pi=rep(1/16,16),
                  states=c(paste0(c("A"),c("A","C","G","T")),paste0(c("C"),c("A","C","G","T")),paste0(c("G"),c("A","C","G","T")),paste0(c("T"),c("A","C","G","T"))))
class(bigger.jc) <- "substitution.model"

bigger.gtr <- list(Q=bigger.gtr,
                   pi=independent.df,
                   states=c(paste0(c("A"),c("A","C","G","T")),paste0(c("C"),c("A","C","G","T")),paste0(c("G"),c("A","C","G","T")),paste0(c("T"),c("A","C","G","T"))))
class(bigger.gtr) <- "substitution.model"

d0.0 <- list(Q=d0.0,
             pi=df,
             states=c(paste0(c("A"),c("A","C","G","T")),paste0(c("C"),c("A","C","G","T")),paste0(c("G"),c("A","C","G","T")),paste0(c("T"),c("A","C","G","T"))))
class(d0.0) <- "substitution.model"

d0.5 <- list(Q=d0.5,
             pi=df,
             states=c(paste0(c("A"),c("A","C","G","T")),paste0(c("C"),c("A","C","G","T")),paste0(c("G"),c("A","C","G","T")),paste0(c("T"),c("A","C","G","T"))))
class(d0.5) <- "substitution.model"

d2.0 <- list(Q=d2.0,
             pi=df,
             states=c(paste0(c("A"),c("A","C","G","T")),paste0(c("C"),c("A","C","G","T")),paste0(c("G"),c("A","C","G","T")),paste0(c("T"),c("A","C","G","T"))))
class(d2.0) <- "substitution.model"

d8.0 <- list(Q=d8.0,
             pi=df,
             states=c(paste0(c("A"),c("A","C","G","T")),paste0(c("C"),c("A","C","G","T")),paste0(c("G"),c("A","C","G","T")),paste0(c("T"),c("A","C","G","T"))))
class(d8.0) <- "substitution.model"

d1000.0 <- list(Q=d1000.0,
                pi=df,
                states=c(paste0(c("A"),c("A","C","G","T")),paste0(c("C"),c("A","C","G","T")),paste0(c("G"),c("A","C","G","T")),paste0(c("T"),c("A","C","G","T"))))
class(d1000.0) <- "substitution.model"

# count substitutions over tree in doublet land
# divide by 2 because rate matrices are normalized to sites but are on site-pairs
count.bigger.jc <- countPriorSubs(bigger.jc,tree)/2

count.bigger.gtr <- countPriorSubs(bigger.gtr,tree)/2

count.d.0.0 <- countPriorSubs(d0.0,tree)/2

count.d.0.5 <- countPriorSubs(d0.5,tree)/2

count.d.2.0 <- countPriorSubs(d2.0,tree)/2

count.d.8.0 <- countPriorSubs(d8.0,tree)/2

count.d.1000.0 <- countPriorSubs(d1000.0,tree)/2

counts <- rbind(count.bigger.jc,count.bigger.gtr,count.d.0.0,count.d.0.5,count.d.2.0,count.d.8.0,count.d.1000.0)

dist(counts)

# comparing in doublet land is hard because things are normalized per site, what if we compared in site land?
site.count.bigger.jc <- doubletCounts2siteCounts(count.bigger.jc)

site.count.bigger.gtr <- doubletCounts2siteCounts(count.bigger.gtr)

site.count.d.0.0 <- doubletCounts2siteCounts(count.d.0.0)

site.count.d.0.5 <- doubletCounts2siteCounts(count.d.0.5)

site.count.d.2.0 <- doubletCounts2siteCounts(count.d.2.0)

site.count.d.8.0 <- doubletCounts2siteCounts(count.d.8.0)

site.count.d.1000.0 <- doubletCounts2siteCounts(count.d.1000.0)

site.counts <- rbind(site.count.bigger.jc,site.count.bigger.gtr,site.count.d.0.0,site.count.d.0.5,site.count.d.2.0,site.count.d.8.0,site.count.d.1000.0)

dist(site.counts)

pca <- prcomp(site.counts)

pca.coords <- site.counts %*% pca$rotation

plot(NULL,NULL,xlim=range(pca.coords[,1])*c(1,1.05),ylim=range(pca.coords[,2])*c(1.05,1.05),xlab="PC1",ylab="PC2")
text(c("JC","GTR",paste0("d=",c(0,0.5,2,8,1000))),x=pca.coords[,1],y=pca.coords[,2])

plot(NULL,NULL,xlim=range(pca.coords[,1])*c(1,1.05),ylim=range(pca.coords[,2])*c(1.05,1.05),xlab="PC1",ylab="PC2")
points(pca.coords[,1:2],col=c("red","black",viridis_pal()(5)),pch=16,cex=2)
text(c("JC","GTR",paste0("d=",c(0,0.5,2,8,1000))),x=pca.coords[,1],y=pca.coords[,2])

