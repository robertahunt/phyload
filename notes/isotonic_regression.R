library(viridis)

set.seed(42)

nboot <- 1000

colors <- viridis_pal()

fit.r <- function(y,n.i,n.e,y.is.decreasing) {
  # recover()
  # function requires monotonically increasing y, so if y is decreasing, we model -y instead
  if ( y.is.decreasing) {
    y <- -y
  }
  fn <- function(r) {
    n_eff <- n.i + r * n.e
    # key <- order(n_eff)
    # n_eff <- n_eff[key]
    # y <- y[key]
    fit <- isoreg(x=n_eff,y=y)
    sum(residuals(fit)^2)
  }
  
  # opts <- list( "algorithm" = "NLOPT_LN_PRAXIS",
  #               "xtol_rel"  = 1.0e-6,
  #               "maxeval"   = 1000 )
  # nloptr(x0=0.5,eval_f=fn,opts=opts)$solution
  
  optimise(fn,c(-1,1))$minimum
}

all.d <- c(0,0.5,2,8,1000)

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

pdf(file="notes/figures/relative_worth.pdf",width=4,height=6)
  par(xpd=TRUE,mfrow=c(2,1),mai=c(0.5,0.5,0.01,0.01),omi=c(0.5,0.5,0.01,0.01))
  historidge(bootstrapped.accuracy,colors=viridis_pal()(5),ylab="accuracy",xlab="")
  legend("bottomleft",legend=rev(paste0("d = ",c(0,0.5,2,8,1000))),fill=rev(viridis_pal()(5)),bty="n",border=NA)
  historidge(bootstrapped.precision,colors=viridis_pal()(5),xlab="relative worth",ylab="precision")
dev.off()

pdf(file="notes/figures/relative_worth.pdf",width=7,height=4)
  par(xpd=TRUE,mfrow=c(1,2),mai=c(0.01,0.01,0.01,0.01),omi=c(0.75,0.75,0.1,0.01))
  boxplot(bootstrapped.accuracy,names=paste0("d=",all.d),col=viridis_pal()(5),outline=FALSE,ylim=c(0.35,1))
  mtext("relative worth",side=2,line=2.5)
  mtext("accuracy-based",side=1,line=2.5)
  boxplot(bootstrapped.precision,names=paste0("d=",all.d),col=viridis_pal()(5),yaxt="n",outline=FALSE,ylim=c(0.35,1))
  mtext("precision-based",side=1,line=2.5)
dev.off()

colors <- viridis_pal()(5)
plot.d <- c(1/1000,0.5,2,8,1000)

plot(NULL,NULL,xlim=c(0.001,1000),ylim=c(0,1),log="x")
for (i in 1:5) {
  lines(x=rep(plot.d[i],2),y=quantile(bootstrapped.accuracy[[i]],c(0.025,0.975)),lwd=2,col="red")
}
lines(plot.d,mle.accuracy,col="red")

for (i in 1:5) {
  lines(x=rep(plot.d[i],2),y=quantile(bootstrapped.precision[[i]],c(0.025,0.975)),lwd=2,col="blue")
}
lines(plot.d,mle.precision,col="blue")
