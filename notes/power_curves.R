library(viridis)

# for plotting
all.d <- c(0.0,0.5,2.0,8.0,1000.0)

# False positives GY93

pvals <- read.csv("~/Downloads/2020-01-06_csv/agg_G93_pvalue.csv",stringsAsFactors=FALSE)

true.pos <- sapply(all.d, function(d){
  true.pos.p <- pvals$value[pvals$X.d. == d & pvals$X.n_e. > 0]
  return(sum(true.pos.p < 0.05)/length(true.pos.p))
  # false.pos.p <- pvals$value[pvals$X.d. == d & pvals$X.n_e. == 0]
  # c(sum(true.pos.p < 0.05)/length(true.pos.p),sum(false.pos.p < 0.05)/length(false.pos.p))
})

null.dist <- pvals$value[pvals$X.n_e. == 0 | pvals$X.d. == 0]
hist(null.dist,breaks=seq(0,1,0.05),xlab="p-value (MI max)",main="n_epi = 0, all d")
sum(null.dist < 0.05)/length(null.dist)

# False positives kurt

pvals <- read.csv("~/Downloads/2020-01-06_csv/agg_kurtosis_pvalue.csv",stringsAsFactors=FALSE)

true.pos <- sapply(all.d, function(d){
  true.pos.p <- pvals$value[pvals$X.d. == d & pvals$X.n_e. > 0]
  return(sum(true.pos.p < 0.05)/length(true.pos.p))
  # false.pos.p <- pvals$value[pvals$X.d. == d & pvals$X.n_e. == 0]
  # c(sum(true.pos.p < 0.05)/length(true.pos.p),sum(false.pos.p < 0.05)/length(false.pos.p))
})

null.dist <- pvals$value[pvals$X.n_e. == 0 | pvals$X.d. == 0]
hist(null.dist,breaks=seq(0,1,0.05),xlab="p-value (MI max)",main="n_epi = 0, all d")
sum(null.dist < 0.05)/length(null.dist)

# False positives max

pvals <- read.csv("~/Downloads/2020-01-06_csv/agg_max_pvalue.csv",stringsAsFactors=FALSE)

true.pos <- sapply(all.d, function(d){
  true.pos.p <- pvals$value[pvals$X.d. == d & pvals$X.n_e. > 0]
  return(sum(true.pos.p < 0.05)/length(true.pos.p))
  # false.pos.p <- pvals$value[pvals$X.d. == d & pvals$X.n_e. == 0]
  # c(sum(true.pos.p < 0.05)/length(true.pos.p),sum(false.pos.p < 0.05)/length(false.pos.p))
})

null.dist <- pvals$value[pvals$X.n_e. == 0 | pvals$X.d. == 0]
hist(null.dist,breaks=seq(0,1,0.05),xlab="p-value (MI max)",main="n_epi = 0, all d")
sum(null.dist < 0.05)/length(null.dist)


## Plot POWER CURVES!!!

pdf("notes/figures/power_curves.pdf",width=6,height=2.5)

  par(mfrow=c(1,3),mai=c(0.7,0.55,0.25,0.01),omi=rep(0.01,4))
  windows.x <- seq(0,1,0.1)
  x.avg <- (windows.x[-1] + windows.x[-length(windows.x)])/2
  cols <- viridis_pal()(5)
  
  x <- sort(c(0,1,rep(seq(0.1,0.9,0.1),2)))
  
  # GY93
  cat("GY93 max power:\n")
  plot(NULL,NULL,xlim=c(-0.1,1),ylim=c(0,1),xlab="proportion of epistatic sites",ylab="power",main="GY93")
  for (i in 1:5) {
    d <- all.d[i]
    
    pvals <- read.csv("~/Downloads/2020-01-06_csv/agg_G93_pvalue.csv",stringsAsFactors=FALSE)
    pvals <- pvals[pvals$X.d. == d,]
    detect <- pvals$value < 0.05
    prop_epi <- pvals$X.n_e./(pvals$X.n_e.+pvals$X.n_i.)
    
    windows.y <- sapply(2:length(windows.x),function(i){
      xl <- windows.x[i-1]
      xh <- ifelse(i == length(windows.x),1.1,windows.x[i])
      indices <- prop_epi < xh & prop_epi >= xl
      sum(detect[indices])/sum(indices)
    })
    cat("  ",max(windows.y),"\n")
    # lines(x.avg,windows.y,col=cols[i],lwd=2)
    
    y <- as.numeric(rbind(windows.y,windows.y))
      
    lines(x,y,col=cols[i],lwd=2)
    
  }

  # kurt
  cat("kurt max power:\n")
  plot(NULL,NULL,xlim=c(-0.1,1),ylim=c(0,1),xlab="proportion of epistatic sites",ylab="power",main="kurtosis(MI)")
  for (i in 1:5) {
    d <- all.d[i]
    
    pvals <- read.csv("~/Downloads/2020-01-06_csv/agg_kurtosis_pvalue.csv",stringsAsFactors=FALSE)
    pvals <- pvals[pvals$X.d. == d,]
    detect <- pvals$value < 0.05
    prop_epi <- pvals$X.n_e./(pvals$X.n_e.+pvals$X.n_i.)
    
    windows.y <- sapply(2:length(windows.x),function(i){
      xl <- windows.x[i-1]
      xh <- ifelse(i == length(windows.x),1.1,windows.x[i])
      indices <- prop_epi < xh & prop_epi >= xl
      sum(detect[indices])/sum(indices)
    })
    cat("  ",max(windows.y),"\n")
    
    # lines(x.avg,windows.y,col=cols[i],lwd=2)
    
    y <- as.numeric(rbind(windows.y,windows.y))
    
    lines(x,y,col=cols[i],lwd=2)
    
  }
  legend("topleft",legend=rev(paste0("d = ",all.d)),fill=rev(cols),bty="n",border=NA,cex=1)
  
  # max
  cat("max max power:\n")
  plot(NULL,NULL,xlim=c(-0.1,1),ylim=c(0,1),xlab="proportion of epistatic sites",ylab="power",main="max(MI)")
  for (i in 1:5) {
    d <- all.d[i]
    
    pvals <- read.csv("~/Downloads/2020-01-06_csv/agg_max_pvalue.csv",stringsAsFactors=FALSE)
    pvals <- pvals[pvals$X.d. == d,]
    detect <- pvals$value < 0.05
    prop_epi <- pvals$X.n_e./(pvals$X.n_e.+pvals$X.n_i.)
    
    windows.y <- sapply(2:length(windows.x),function(i){
      xl <- windows.x[i-1]
      xh <- ifelse(i == length(windows.x),1.1,windows.x[i])
      indices <- prop_epi < xh & prop_epi >= xl
      sum(detect[indices])/sum(indices)
    })
    cat("  ",max(windows.y),"\n")
    
    # lines(x.avg,windows.y,col=cols[i],lwd=2)
    
    y <- as.numeric(rbind(windows.y,windows.y))
    
    lines(x,y,col=cols[i],lwd=2)
    
  }
dev.off()









# 
# 
# pvals <- read.csv("~/Downloads/2020-01-06_csv/agg_max_pvalue.csv",stringsAsFactors=FALSE)
# dists <- read.csv("~/Downloads/2020-01-06_csv/agg_rf.mean.csv",stringsAsFactors=FALSE)
# 
# pdf("~/Downloads/pvalue_and_rf.pdf",width=6,height=10)
# par(mfrow=c(5,3),mai=c(0.5,0.5,0.4,0.05))
# for (d in all.d) {
#   pvals <- read.csv("~/Downloads/2020-01-06_csv/agg_max_pvalue.csv",stringsAsFactors=FALSE)
#   pvals <- pvals[pvals$X.d. == d,]
#   dists <- read.csv("~/Downloads/2020-01-06_csv/agg_rf.mean.csv",stringsAsFactors=FALSE)
#   dists <- dists[dists$X.d. == d,]
#   
#   # p vs prop
#   plot(pvals$X.n_e./(pvals$X.n_e.+pvals$X.n_i.),-log(pvals$value,base=10),xlab="proportion of epistatic sites",ylab="-log10(p)",main=paste0("d = ",d),pch=16,col="#66666690")
#   fit <- lowess(-log(pvals$value,base=10) ~ pvals$X.n_e./(pvals$X.n_e.+pvals$X.n_i.))
#   lines(fit$x,fit$y,col="red",lwd=2)
#   abline(h=-log(0.05,base=10),lty=2,col="red")
#   
#   # RF vs prop
#   plot(dists$X.n_e./(dists$X.n_e.+dists$X.n_i.),dists$value,xlab="proportion of epistatic sites",ylab="RF distance",main=paste0("d = ",d),pch=16,col="#66666690")
#   fit <- lowess(dists$value ~ dists$X.n_e./(dists$X.n_e.+dists$X.n_i.))
#   lines(fit$x,fit$y,col="red",lwd=2)
#   abline(h=-log(0.05,base=10),lty=2,col="red")
# 
#   # rf vs p
#   plot(dists$value,-log(pvals$value,base=10),xlab="RF dist",ylab="-log10(p)",main=paste0("d = ",d),pch=16,col="#66666690")
#   fit <- lowess(-log(pvals$value,base=10) ~ dists$value)
#   lines(fit$x,fit$y,col="red",lwd=2)
#   abline(h=-log(0.05,base=10),lty=2,col="red")
# }
# dev.off()
# 
# pdf("~/Downloads/pvalue_and_rf.pdf",width=6,height=10)
# par(mfrow=c(5,3),mai=c(0.5,0.5,0.4,0.05))
# for (d in all.d) {
#   pvals <- read.csv("~/Downloads/2020-01-06_csv/agg_max_pvalue.csv",stringsAsFactors=FALSE)
#   pvals <- pvals[pvals$X.d. == d,]
#   pvals <- pvals[pvals$X.n_i. + pvals$X.n_e. == 400,]
#   dists <- read.csv("~/Downloads/2020-01-06_csv/agg_rf.mean.csv",stringsAsFactors=FALSE)
#   dists <- dists[dists$X.d. == d,]
#   dists <- dists[dists$X.n_i. + dists$X.n_e. == 400,]
#   
#   # p vs prop
#   plot(pvals$X.n_e./(pvals$X.n_e.+pvals$X.n_i.),-log(pvals$value,base=10),xlab="proportion of epistatic sites",ylab="-log10(p)",main=paste0("d = ",d),pch=16,col="#66666690")
#   fit <- lowess(-log(pvals$value,base=10) ~ pvals$X.n_e./(pvals$X.n_e.+pvals$X.n_i.))
#   lines(fit$x,fit$y,col="red",lwd=2)
#   abline(h=-log(0.05,base=10),lty=2,col="red")
#   
#   # RF vs prop
#   plot(dists$X.n_e./(dists$X.n_e.+dists$X.n_i.),dists$value,xlab="proportion of epistatic sites",ylab="RF distance",main=paste0("d = ",d),pch=16,col="#66666690")
#   fit <- lowess(dists$value ~ dists$X.n_e./(dists$X.n_e.+dists$X.n_i.))
#   lines(fit$x,fit$y,col="red",lwd=2)
#   abline(h=-log(0.05,base=10),lty=2,col="red")
#   
#   # rf vs p
#   plot(dists$value,-log(pvals$value,base=10),xlab="RF dist",ylab="-log10(p)",main=paste0("d = ",d),pch=16,col="#66666690")
#   fit <- lowess(-log(pvals$value,base=10) ~ dists$value)
#   lines(fit$x,fit$y,col="red",lwd=2)
#   abline(h=-log(0.05,base=10),lty=2,col="red")
# }
# dev.off()
