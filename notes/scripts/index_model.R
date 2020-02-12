library(np)
library(viridis)

# Point estimates

error <- read.csv("~/Downloads/2020-01-06_csv/agg_rf.mean.csv",stringsAsFactors=FALSE)
all.fits.error <- lapply(c(0.0,0.5,2.0,8.0,1000.0),function(d){
  D <- error[error$X.d. == d,]
  
  X <- cbind(rep(1,dim(D)[1]),D$X.n_i.,D$X.n_e.)
  y <- log(D$value)
  
  bw <- npindexbw(xdat=X,ydat=y)
  npindex(bws=bw,gradients=TRUE)
  
})

precision <- read.csv("~/Downloads/2020-01-06_csv/agg_mrc_percent_resolved.csv",stringsAsFactors=FALSE)
precision$value[precision$value > (1 - 1e-8)] <- (1 - 1e-8)
all.fits.precision <- lapply(c(0.0,0.5,2.0,8.0,1000.0),function(d){
  D <- precision[precision$X.d. == d,]
  
  X <- cbind(rep(1,dim(D)[1]),D$X.n_i.,D$X.n_e.)
  y <- LaplacesDemon::logit(1 - D$value)
  
  bw <- npindexbw(xdat=X,ydat=y)
  npindex(bws=bw,gradients=TRUE)
  
})

error.worth <- unlist(lapply(all.fits.error,function(fit){
  beta <- coef(fit)
  beta[3]/beta[2]
}))

error.worth.ci <- do.call(rbind,lapply(all.fits.error,function(fit){
  samples <- LaplacesDemon::rmvn(1e6,fit$beta[-1],fit$betavcov[-1,-1])
  quantile(samples[,2]/samples[,1],prob=c(0.025,0.975))
}))

precision.worth <- unlist(lapply(all.fits.precision,function(fit){
  beta <- coef(fit)
  beta[3]/beta[2]
}))

all.d <- c(0.0,0.5,2.0,8.0,1000.0)

plot(all.d,error.worth,type="l")
lines(all.d,precision.worth,type="l")


# Jack-knifing
nrep <- 100
nsubsample <- 300

all.fits.error.jack <- lapply(c(0.0,0.5,2.0,8.0,1000.0),function(d){
  D <- error[error$X.d. == d,]
  X <- cbind(rep(1,dim(D)[1]),D$X.n_i.,D$X.n_e.)
  y <- LaplacesDemon::logit(1 - D$value)
  
  sapply(1:nrep,function(i){
    indices <- sample(1:dim(D)[1],300)
    X_rep <- X[indices,]
    y_rep <- y[indices]
    
    bw <- npindexbw(xdat=X_rep,ydat=y_rep)
    fit <- npindex(bws=bw)
    fit$beta[3]/fit$beta[2]
  })
  
})

names(all.fits.error.jack) <- all.d
boxplot(all.fits.error.jack,xlab="d",ylab="relative worth")

all.fits.precision.jack <- lapply(c(0.0,0.5,2.0,8.0,1000.0),function(d){
  D <- precision[precision$X.d. == d,]
  X <- cbind(rep(1,dim(D)[1]),D$X.n_i.,D$X.n_e.)
  y <- log(D$value)
  
  sapply(1:nrep,function(i){
    indices <- sample(1:dim(D)[1],300)
    X_rep <- X[indices,]
    y_rep <- y[indices]
    
    bw <- npindexbw(xdat=X_rep,ydat=y_rep)
    fit <- npindex(bws=bw)
    fit$beta[3]/fit$beta[2]
  })
  
})

names(all.fits.precision.jack) <- all.d
boxplot(all.fits.precision.jack,xlab="d",ylab="relative worth")



comparogram(all.fits.error.jack,all.fits.precision.jack,left.color=myorange,right.color=myblue,xlab="d",names=all.d,ylab="relative worth")
legend("topright",fill=c(myorange,myblue),legend=c("error","precision"),bty="n",border=NA)
