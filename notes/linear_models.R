library(MASS)
library(viridis)

applyByD <- function(data,fn,d=c(0.0,0.5,2.0,8.0,1000.0)) {
  x <- sapply(d,function(d_){fn(data$value[data$X.d. == d_])})
  names(x) <- d
  return(x)
}

toMatrices <- function(data,d=c(0.0,0.5,2.0,8.0,1000.0)) {
  x <- lapply(d,function(d_){
    tmp <- matrix(nrow=26,ncol=26)
    for (n in 1:675) {
      i <- data$X.n_i.[data$X.d. == d_][n]/16 + 1
      j <- data$X.n_e.[data$X.d. == d_][n]/16 + 1
      tmp[i,j] <- data$value[data$X.d. == d_][n]
    }
    rownames(tmp) <- seq(0,400,16)
    colnames(tmp) <- seq(0,400,16)
    return(tmp)
  }) 
  names(x) <- d
  return(x)
}

plotMainDiagonal <- function(dmat) {
  y <- diag(t(apply(dmat,2,rev)))
  niid <- seq(400,0,-16)
  nepi <- seq(0,400,16)
  x <- nepi/(nepi + niid)
  plot(x,y,xlab="prop epi",ylab="value",main="")
}

pointsMainDiagonal <- function(dmat,col,pch=16,regress=TRUE) {
  y <- diag(t(apply(dmat,2,rev)))
  niid <- seq(400,0,-16)
  nepi <- seq(0,400,16)
  x <- nepi/(nepi + niid)
  points(x,y,pch=pch,col=col)
  if (regress) {
    abline(lm(y~x),lty=2,col=col)
  }
}

getMainDiagonal <- function(dmat) {
  y <- diag(t(apply(dmat,2,rev)))
  niid <- seq(400,0,-16)
  nepi <- seq(0,400,16)
  x <- nepi/(nepi + niid)
  return(cbind(x,y))
}

plotAllMainDiagonals <- function(alldmat,ylab,cols,legend.location,legend.cex,verbose=FALSE) {
  xy <- lapply(1:length(alldmat),function(i){
    getMainDiagonal(alldmat[[i]])
  })
  y <- lapply(xy,function(xy_){xy_[,2]})
  plot(NULL,NULL,xlab="proportion epistatic",ylab=ylab,main="",xlim=c(0,1),ylim=range(y))
  for (i in 1:length(xy)) {
    points(xy[[i]][,"x"],xy[[i]][,"y"],col=cols[i],pch=16)
    abline(lm(xy[[i]][,"y"] ~ xy[[i]][,"x"]),col=cols[i],lwd=2)
    if (verbose) {
      print(summary(lm(xy[[i]][,"y"] ~ xy[[i]][,"x"])))
    }
  }
  legend(legend.location,legend=names(alldmat),fill=cols,bty="n",border=NA,cex=legend.cex)
}

extractAllContrasts <- function(alldmat,fxn=c("+","-","*","/")) {
  n <- dim(alldmat[[1]])[1]
  lapply(alldmat,function(dm){
    mat <- matrix(NA,n-1,n-1)
    for (i in 1:(n-1)) {
      for (j in 1:(n-1)) {
        if ( fxn == "+" ) {
          mat[i,j] <- (dm[i,j+1] - dm[i,j]) + (dm[i+1,j] - dm[i,j])
        } else if ( fxn == "-" ) {
          mat[i,j] <- (dm[i,j+1] - dm[i,j]) - (dm[i+1,j] - dm[i,j])
        } else if ( fxn == "*" ) {
          mat[i,j] <- (dm[i,j+1] - dm[i,j]) * (dm[i+1,j] - dm[i,j])
        } else if ( fxn == "/" ) {
          mat[i,j] <- (dm[i,j+1] - dm[i,j]) / (dm[i+1,j] - dm[i,j])
        } else {
          stop("Invalid function")
        }
        
      }
    }
    return(mat)
  })
}

centralDifferenceDeltaRatios <- function(alldmat) {
  n <- dim(alldmat[[1]])[1]
  lapply(alldmat,function(dm){
    mat <- matrix(NA,n-2,n-2)
    for (i in 2:(n-1)) {
      for (j in 2:(n-1)) {
        mat[i-1,j-1] <- (dm[i,j+1] - dm[i,j-1]) / (dm[i+1,j] - dm[i-1,j])
      }
    }
    return(as.numeric(mat))
  })
}

extractAllProportionateError <- function(alldmat) {
  n <- dim(alldmat[[1]])[1]
  lapply(alldmat,function(dm){
    mat <- matrix(NA,n-1,n-1)
    for (i in 1:(n-1)) {
      for (j in 1:(n-1)) {
        mat[i,j] <- (((dm[i,j+1] - dm[i,j]) - (dm[i+1,j] - dm[i,j])))/(dm[i+1,j] - dm[i,j])
      }
    }
    return(mat)
  })
}

getAllForwardDifferences <- function(alldmat,direction=c("right/left","up/down")) {
  n <- dim(alldmat[[1]])[1]
  lapply(alldmat,function(dm){
    if ( direction == "right" || direction == "left" ) {
      mat <- sapply(1:(n-1),function(j){
        dm[,j+1] - dm[,j]
      })
    } else if ( direction == "up" || direction == "down" ) {
      mat <- sapply(1:(n-1),function(i){
        dm[i+1,] - dm[i,]
      })
    } else {
      stop("Invalid direction")
    }
    return(t(mat))
  })
}

baseline <- function(alldmat) {
  do.call(rbind,lapply(alldmat,function(dm){
    dm[,1]
  }))
}

data <- read.csv("~/Downloads/2020-01-06_csv/agg_rf.mean.csv",stringsAsFactors=FALSE)

applyByD(data,mean)
applyByD(data,median)

dm <- toMatrices(data)

boxplot(lapply(1:5,function(i){baseline(dm)[i,]}))




delta.e <- getAllForwardDifferences(dm,"right")
delta.i <- getAllForwardDifferences(dm,"down")

nsites <- t(sapply(1:25,function(i){rep(8 + 16*(i-1),26)}))
nsites <- as.numeric(nsites)

groups <- as.factor(as.numeric(sapply(1:26,function(i){rep(i,25)})))

worth.error1 <- sapply(1:5,function(index){
  di <- as.numeric(delta.i[[index]])
  de <- as.numeric(delta.e[[index]])
  fit.i <- rlm(di ~ nsites + groups -1)
  fit.e <- rlm(de ~ nsites + groups -1)
  return(coef(fit.e)["nsites"]/coef(fit.i)["nsites"])
})

worth.error2 <- sapply(1:5,function(index){
  di <- as.numeric(delta.i[[index]])
  de <- as.numeric(delta.e[[index]])
  fit.i <- rlm(di ~ nsites)
  fit.e <- rlm(de ~ nsites)
  return(coef(fit.e)["nsites"]/coef(fit.i)["nsites"])
})

plot(nsites,delta.i[[2]])
abline(rlm(as.numeric(delta.i[[2]]) ~ as.numeric(nsites)))

summary(rlm(as.numeric(delta.i[[2]]) ~ as.numeric(nsites)))

summary()
summary()


plotAllMainDiagonals(dm,"p-val",viridis_pal(option="viridis")(5),"topleft",0.5,verbose=FALSE)



n.iid   <- data$X.n_i.
n.epi   <- data$X.n_e.
d       <- data$X.d.
d_logit <- LaplacesDemon::invlogit(d)
badness <- data$value

real.scale <- lm(badness ~ n.iid + n.epi + d_logit:n.epi)
real.scale$coefficients["n.epi"]/real.scale$coefficients["n.iid"]
real.scale$coefficients["n.epi:d_logit"]/real.scale$coefficients["n.iid"]








real.scale <- lm(badness ~ n.iid + n.epi + d:n.epi)
real.scale$coefficients["n.epi"]/real.scale$coefficients["n.iid"]
real.scale$coefficients["n.epi:d"]/real.scale$coefficients["n.iid"]

real.scale.robust <- rlm(badness ~ n.iid + n.epi + d:n.epi)
real.scale.robust$coefficients["n.epi"]/real.scale.robust$coefficients["n.iid"]
real.scale.robust$coefficients["n.epi:d"]/real.scale.robust$coefficients["n.iid"]


n.iid   <- data$X.n_i.
n.epi   <- data$X.n_e.
d       <- data$X.d.
badness <- log(data$value)

log.scale.robust <- rlm(badness ~ n.iid + n.epi + d:n.epi)
log.scale.robust$coefficients["n.epi"]/log.scale.robust$coefficients["n.iid"]
log.scale.robust$coefficients["n.epi:d"]/log.scale.robust$coefficients["n.iid"]

coef.epi <- log.scale.robust$coefficients["n.epi"]/log.scale.robust$coefficients["n.iid"]
coef.d   <- log.scale.robust$coefficients["n.epi:d"]/log.scale.robust$coefficients["n.iid"]


effectiveSequenceLength <- function(n.iid,n.epi,d,coef.epi,coef.d) {
  n.iid + n.epi * coef.epi + n.epi * coef.epi * coef.d * d
}

effectiveSequenceLength(200,200,0.5,coef.epi,coef.d)
effectiveSequenceLength(200,200,2,coef.epi,coef.d)
effectiveSequenceLength(200,200,8,coef.epi,coef.d)







data <- read.csv("~/Downloads/tree_badness/agg_rf.mean.csv")


n.iid   <- data$X.n_i.
n.epi   <- data$X.n_e.
p.doublet <- data$X.d.
p.doublet[data$X.d. == 0.5] <- 0.03191064
p.doublet[data$X.d. == 2.0] <- 0.1164907
p.doublet[data$X.d. == 8.0] <- 0.3452926
n.epi.adj <- n.epi * 0.5 * (2 - p.doublet)

badness <- log(data$value)

mechanistic <- rlm(badness ~ n.iid + n.epi.adj)
mechanistic$coefficients["n.epi.adj"]/mechanistic$coefficients["n.iid"]

summary(mechanistic)



n.iid   <- data$X.n_i.
n.epi   <- data$X.n_e.
d       <- data$X.d.
n.epi.d.05 <- n.epi * (d == 0.5)
n.epi.d.2 <- n.epi * (d == 2.0)
n.epi.d.8 <- n.epi * (d == 8.0)
badness <- log(data$value)

per.d <- rlm(badness ~ n.iid + n.epi.d.05 + n.epi.d.2 + n.epi.d.8)
per.d$coefficients["n.epi.d.05"]/per.d$coefficients["n.iid"]
per.d$coefficients["n.epi.d.2"]/per.d$coefficients["n.iid"]
per.d$coefficients["n.epi.d.8"]/per.d$coefficients["n.iid"]

mechanistic$coefficients["n.epi.adj"]/mechanistic$coefficients["n.iid"] * 0.5 * (2 - 0.03191064)
mechanistic$coefficients["n.epi.adj"]/mechanistic$coefficients["n.iid"] * 0.5 * (2 - 0.1164907)
mechanistic$coefficients["n.epi.adj"]/mechanistic$coefficients["n.iid"] * 0.5 * (2 - 0.3452926)











library(viridis)



rf.05 <- read.csv("~/Downloads/tree_badness/agg_rf.q05.csv")
rf.95 <- read.csv("~/Downloads/tree_badness/agg_rf.q95.csv")

rf.05 <- rf.05[rf.05$X.d. == 2,]
rf.95 <- rf.95[rf.95$X.d. == 2,]

n.iid    <- rf.05$X.n_i.
n.epi    <- rf.05$X.n_e.
rf.width <- rf.95$value - rf.05$value

magma <- viridis_pal(option="magma")(100)

heats <- cut(rf.width,breaks=100,include.lowest = TRUE,labels = FALSE)
heats <- viridis_pal(option="magma")(100)[heats]

nsites <- seq(4,400,16)

par(mfrow=c(1,2))

plot(NULL,NULL,xlim=c(0,26),ylim=c(0,26),xlab="n iid",ylab="n epi",xaxt="n",yaxt="n")
axis(1,at=seq(1,25,4),labels=nsites[seq(1,25,4)])
axis(2,at=seq(1,25,4),labels=nsites[seq(1,25,4)])
for (i in 1:25) {
  for (j in 1:25) {
    
    xmin <- i - 1
    xmax <- i
    ymin <- j - 1
    ymax <- j
    
    key <- which(rf.05$X.n_i. == nsites[i] & rf.05$X.n_e. == nsites[j])
    
    polygon(x=c(xmin,xmin,xmax,xmax),y=c(ymin,ymax,ymax,ymin),border=NA,col=heats[key])
    
  }
}

plot(NULL,NULL,xlim=c(0,1),ylim=c(0,100),xlab="",ylab="",xaxt="n",bty="n",yaxt="n")
rf <- quantile(rf.width,(seq(0,100,20))/100)
axis(2,at=seq(0,100,20),labels=rf)

xmin <- 0
xmax <- 1
for (i in 1:length(heats)) {
  ymin <- i - 1
  ymax <- i
  
  polygon(x=c(xmin,xmin,xmax,xmax),y=c(ymin,ymax,ymax,ymin),border=NA,col=magma[i])
  
}

summary(rlm(rf.width ~ n.iid + n.epi))















rf.05 <- read.csv("~/Downloads/tree_badness/agg_rf.q05.csv")
rf.95 <- read.csv("~/Downloads/tree_badness/agg_rf.q95.csv")

n.iid    <- rf.05$X.n_i.
n.epi    <- rf.05$X.n_e.
d        <- rf.05$X.d.
rf.width <- rf.95$value - rf.05$value

summary(rlm(rf.width ~ n.iid + n.epi + d:n.epi))
