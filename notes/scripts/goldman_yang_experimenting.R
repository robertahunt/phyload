library(LaplacesDemon)

JSD <- function(P,Q) {
  M <- (P + Q)/2
  0.5*(sum(P * log(P/M,base=2))) + 0.5*(sum(Q * log(Q/M,base=2)))
}

set.seed(47)

er <- c(2,8,1,1.5,7,1)
bf1 <- c(0.4,0.2,0.3,0.1)

phy <- rtree(50,FALSE,br=rexp)
phy$edge.length <- phy$edge.length/sum(phy$edge.length) * 0.5


aln1 <- lapply(1:100,function(i){simSeq(phy,Q=er,bf=bf1,l=800)})
gy1 <- unlist(lapply(aln1,unconstrainedLikelihood))

bf <- rdirichlet(10000,c(1,1,1,1))
divs <- sapply(1:10000,function(i){JSD(bf[i,],bf1)})
bf <- bf[divs > 0.02 & divs < 0.023,]

gy_other <- vector("list",dim(bf)[1])
divs <- divs[divs > 0.02 & divs < 0.023]

for (i in 1:dim(bf)[1]) {
  bf2 <- bf[i,]
  
  aln2 <- lapply(1:100,function(i){simSeq(phy,Q=er,bf=bf2,l=800)})
  gy_other[[i]] <- unlist(lapply(aln2,unconstrainedLikelihood))
  
}

gy_dist <- unlist(lapply(gy_other,mean))

hist(gy_dist,breaks=20,col="#66666690",border=NA)
abline(v=mean(gy1),col="red")

discernable <- sapply(1:length(gy_other),function(i){
  (sum(gy1 < quantile(gy_other[[i]],0.025)) + sum(gy1 > quantile(gy_other[[i]],0.975)))/100
})

hist(discernable,breaks=20,col="#66666690",border=NA,main="",xlab="Pr(detect difference)")

plot(divs,discernable)
