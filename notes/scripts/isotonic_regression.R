library(viridis)

set.seed(42)

nboot <- 1000

colors <- viridis_pal()(5)

all.d <- c(0,0.5,2,8,1000)

fit.r <- function(y,n.i,n.e,y.is.decreasing) {
  # recover()
  # function requires monotonically increasing y, so if y is decreasing, we model -y instead
  if ( y.is.decreasing) {
    y <- -y
  }
  fn <- function(r) {
    n_eff <- n.i + r * n.e
    fit <- isoreg(x=n_eff,y=y)
    sum(residuals(fit)^2)
  }
  # ln_n <- log(length(y))
  # fn_reg <- function(r) {
  #   n_eff <- n.i + r * n.e
  #   fit <- isoreg(x=n_eff,y=y)
  #   ln_n * length(fit$iKnots) - 2*sum(dnorm(residuals(fit),log=TRUE))
  # }
  optimise(fn,c(-1,1))$minimum
}

# fit.r <- function(y,n.i,n.e,y.is.decreasing) {
#   recover()
#   # function requires monotonically increasing y, so if y is decreasing, we model -y instead
#   if ( y.is.decreasing) {
#     y <- -y
#   }
#   fn <- function(r) {
#     n_eff <- n.i + r * n.e
#     # knots <- quantile(n_eff,seq(0.2,0.8,0.2))
#     knots <- seq(min(n_eff),max(n_eff),length.out=21)[-c(1,21)]
#     X <- iSpline(n_eff,knots=knots,degree=1)
#     fit <- rlm(y ~ X)
#     sum(abs(residuals(fit)))
#   }
#   optimise(fn,c(-1,1))$minimum
# }

mle.accuracy <- sapply(all.d,function(d){
  data <- read.csv("~/Downloads/2020-01-06_csv/agg_rf.mean.csv",stringsAsFactors=FALSE)
  data <- data[data$X.d. == d,]
  fit.r(data$value,data$X.n_i.,data$X.n_e.,y.is.decreasing=TRUE)
})

bootstrapped.accuracy <- lapply(all.d,function(d){
  data <- read.csv("~/Downloads/2020-01-06_csv/agg_rf.mean.csv",stringsAsFactors=FALSE)
  data <- data[data$X.d. == d,]
  
  sapply(1:nboot,function(i) {
    idx <- sample.int(dim(data)[1],replace=TRUE)
    data_boot <- data[idx,]
    fit.r(data_boot$value,data_boot$X.n_i.,data_boot$X.n_e.,y.is.decreasing=TRUE)  
  })

})

mle.precision <- sapply(all.d,function(d){
  data <- read.csv("~/Downloads/2020-01-06_csv/agg_mrc_percent_resolved.csv",stringsAsFactors=FALSE)
  data <- data[data$X.d. == d,]
  fit.r(data$value,data$X.n_i.,data$X.n_e.,y.is.decreasing=FALSE)
})

bootstrapped.precision <- lapply(all.d,function(d){
  data <- read.csv("~/Downloads/2020-01-06_csv/agg_mrc_percent_resolved.csv",stringsAsFactors=FALSE)
  data <- data[data$X.d. == d,]
  
  sapply(1:nboot,function(i) {
    idx <- sample.int(dim(data)[1],replace=TRUE)
    data_boot <- data[idx,]
    fit.r(data_boot$value,data_boot$X.n_i.,data_boot$X.n_e.,y.is.decreasing=FALSE)  
  })
  
})

bootstrapped.overconfidence <- lapply(1:5,function(i){
  bootstrapped.precision[[i]]/bootstrapped.accuracy[[i]]
})

pdf(file="notes/figures/relative_worth.pdf",width=7,height=4)
  par(xpd=TRUE,mfrow=c(1,2),mai=c(0.01,0.01,0.01,0.01),omi=c(0.75,0.75,0.1,0.01))
  boxplot(bootstrapped.accuracy,names=paste0("d=",all.d),col=viridis_pal()(5),outline=FALSE,ylim=c(0.35,1))
  mtext("relative worth",side=2,line=2.5)
  mtext("accuracy-based",side=1,line=2.5)
  boxplot(bootstrapped.precision,names=paste0("d=",all.d),col=viridis_pal()(5),yaxt="n",outline=FALSE,ylim=c(0.35,1))
  mtext("precision-based",side=1,line=2.5)
dev.off()

#######
# Sensitivity analyses to convergence standards
#######

# 3 convergence standards:
#  1) anything goes
#  2) discard analyses that don't pass convergence cutoffs
#  3) strict cutoffs
thresh.a <- c(Inf,0.05,0.01)
thresh.p <- c(Inf,1.1,1.01)

sensitivity.accuracy <- lapply(1:3,function(i){
  asdsf <- read.csv("~/Downloads/2020-01-06_csv/agg_asdsf.csv",stringsAsFactors=FALSE)
  psrf <- read.csv("~/Downloads/2020-01-06_csv/agg_psrf.csv",stringsAsFactors=FALSE)
  
  lapply(all.d,function(d){
    data <- read.csv("~/Downloads/2020-01-06_csv/agg_rf.mean.csv",stringsAsFactors=FALSE)
    data <- data[asdsf$value < thresh.a[i] & psrf$value < thresh.p[i],]
    data <- data[data$X.d. == d,]
  
    sapply(1:nboot,function(i) {
      idx <- sample.int(dim(data)[1],replace=TRUE)
      data_boot <- data[idx,]
      fit.r(data_boot$value,data_boot$X.n_i.,data_boot$X.n_e.,y.is.decreasing=TRUE)  
    })
  })
  
})

sensitivity.precision <- lapply(1:3,function(i){
  asdsf <- read.csv("~/Downloads/2020-01-06_csv/agg_asdsf.csv",stringsAsFactors=FALSE)
  psrf <- read.csv("~/Downloads/2020-01-06_csv/agg_psrf.csv",stringsAsFactors=FALSE)
  
  lapply(all.d,function(d){
    data <- read.csv("~/Downloads/2020-01-06_csv/agg_mrc_percent_resolved.csv",stringsAsFactors=FALSE)
    data <- data[asdsf$value < thresh.a[i] & psrf$value < thresh.p[i],]
    data <- data[data$X.d. == d,]
    
    sapply(1:nboot,function(i) {
      idx <- sample.int(dim(data)[1],replace=TRUE)
      data_boot <- data[idx,]
      fit.r(data_boot$value,data_boot$X.n_i.,data_boot$X.n_e.,y.is.decreasing=FALSE)  
    })
  })
  
})

pdf(file="notes/figures/relative_worth_sensitivity.pdf",width=8,height=8)
  # accuracy
  par(xpd=TRUE,mfrow=c(2,3),mai=c(0.01,0.01,0.01,0.01),omi=c(0.75,0.75,0.1,0.01))
  boxplot(sensitivity.accuracy[[1]],names=paste0("d=",all.d),col=viridis_pal()(5),outline=FALSE,ylim=c(0.35,1),xaxt="n")
  mtext("relative worth (accuracy)",side=2,line=2.5)
  boxplot(sensitivity.accuracy[[2]],names=paste0("d=",all.d),col=viridis_pal()(5),outline=FALSE,ylim=c(0.35,1),yaxt="n",xaxt="n")
  boxplot(sensitivity.accuracy[[3]],names=paste0("d=",all.d),col=viridis_pal()(5),outline=FALSE,ylim=c(0.35,1),yaxt="n",xaxt="n")

  # precision
  boxplot(sensitivity.precision[[1]],names=paste0("d=",all.d),col=viridis_pal()(5),outline=FALSE,ylim=c(0.35,1))
  mtext("relative worth (precision)",side=2,line=2.5)
  mtext("no convergence filtering",side=1,line=2.5)
  boxplot(sensitivity.precision[[2]],names=paste0("d=",all.d),col=viridis_pal()(5),outline=FALSE,ylim=c(0.35,1),yaxt="n")
  mtext("convergence filtering",side=1,line=2.5)
  boxplot(sensitivity.precision[[3]],names=paste0("d=",all.d),col=viridis_pal()(5),outline=FALSE,ylim=c(0.35,1),yaxt="n")
  mtext("strict convergence filtering",side=1,line=2.5)
  
dev.off()
