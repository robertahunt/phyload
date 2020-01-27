library(viridis)

applyByD <- function(data,fn,d=c(0.0,0.5,2.0,8.0,1000.0)) {
  x <- sapply(d,function(d_){fn(data$value[data$X.d. == d_])})
  names(x) <- d
  return(x)
}

# For plotting
all.d <- c(0,0.5,2,8,1000)
d.cols <- viridis_pal()(5)
d.cols <- gsub("FF","90",d.cols)
line.cols <- viridis_pal()(5)

# Data
gy93 <- read.csv("~/Downloads/2020-01-06_csv/agg_G93.csv",stringsAsFactors=FALSE)
kurt <- read.csv("~/Downloads/2020-01-06_csv/agg_kurtosis.csv",stringsAsFactors=FALSE)
mmax  <- read.csv("~/Downloads/2020-01-06_csv/agg_max.csv",stringsAsFactors=FALSE)

# sensitity to d
avg.gy93 <- applyByD(gy93,mean)
avg.kurt <- applyByD(kurt,mean)
avg.mmax <- applyByD(mmax,mean)

ratio.gy93 <- sapply(2:5,function(i){avg.gy93[i]/avg.gy93[i-1]})
ratio.kurt <- sapply(2:5,function(i){avg.kurt[i]/avg.kurt[i-1]})
ratio.mmax <- sapply(2:5,function(i){avg.mmax[i]/avg.mmax[i-1]})

mean(ratio.gy93)
mean(ratio.kurt)
mean(ratio.mmax)

# sensitivity to the number of sites
nsites <- gy93$X.n_i. + gy93$X.n_e.
cor(gy93$value,nsites)
cor(kurt$value,nsites)
cor(mmax$value,nsites)

# diagonal-ish values only
gy93 <- gy93[gy93$X.n_i. + gy93$X.n_e. >= 384 & gy93$X.n_i. + gy93$X.n_e. <= 414,]
kurt <- kurt[kurt$X.n_i. + kurt$X.n_e. >= 384 & kurt$X.n_i. + kurt$X.n_e. <= 414,]
mmax <- mmax[mmax$X.n_i. + mmax$X.n_e. >= 384 & mmax$X.n_i. + mmax$X.n_e. <= 414,]

# plot sensitivity to prop epi
pdf("~/git_repos/phyload/notes/figures/test_stat_sensitivity.pdf",width=7,height=2.5)
par(mfrow=c(1,3),mai=c(0.75,0.6,0.05,0.01),omi=rep(0.01,4))
r <- range(gy93$value)
plot(NULL,NULL,xlim=c(0,1),ylim=r,xlab="proportion of epistatic sites",ylab="GY93")
for (i in 1:5) {
  d <- all.d[i]
  points((gy93$X.n_e./(gy93$X.n_e. + gy93$X.n_i.))[gy93$X.d. == d],gy93$value[gy93$X.d. == d],pch=16,col=d.cols[i],cex=1.25)
  fit <- lm(gy93$value[gy93$X.d. == d] ~ (gy93$X.n_e./(gy93$X.n_e. + gy93$X.n_i.))[gy93$X.d. == d])
  abline(fit,lty=1,lwd=3,col=line.cols[i])
}

r <- range(kurt$value)
plot(NULL,NULL,xlim=c(0,1),ylim=r,xlab="proportion of epistatic sites",ylab="kurtosis(MI)")
for (i in 1:5) {
  d <- all.d[i]
  points((kurt$X.n_e./(kurt$X.n_e. + kurt$X.n_i.))[kurt$X.d. == d],kurt$value[kurt$X.d. == d],pch=16,col=d.cols[i],cex=1.25)
  fit <- lm(kurt$value[kurt$X.d. == d] ~ (kurt$X.n_e./(kurt$X.n_e. + kurt$X.n_i.))[kurt$X.d. == d])
  abline(fit,lty=1,lwd=3,col=line.cols[i])
}
legend("topleft",legend=paste0("d = ",all.d),fill=d.cols,bty="n",border=NA)

r <- range(mmax$value)
plot(NULL,NULL,xlim=c(0,1),ylim=r,xlab="proportion of epistatic sites",ylab="max(MI)",log="")
for (i in 1:5) {
  d <- all.d[i]
  points((mmax$X.n_e./(mmax$X.n_e. + mmax$X.n_i.))[mmax$X.d. == d],mmax$value[mmax$X.d. == d],pch=16,col=d.cols[i],cex=1.25)
  fit <- lm(mmax$value[mmax$X.d. == d] ~ (mmax$X.n_e./(mmax$X.n_e. + mmax$X.n_i.))[mmax$X.d. == d])
  abline(fit,lty=1,lwd=3,col=line.cols[i])
}
dev.off()
