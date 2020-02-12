library(viridis)

pvals <- read.csv("~/Downloads/2020-01-06_csv/agg_max_pvalue.csv",stringsAsFactors=FALSE)
rf <- read.csv("~/Downloads/2020-01-06_csv/agg_rf.mean.csv",stringsAsFactors=FALSE)

all.d <- c(0,0.5,2,8,1000)

pdf("notes/figures/power_vs_rf.pdf",width=6,height=2.5)

  par(mfrow=c(1,3),mai=c(0.7,0.55,0.25,0.01),omi=rep(0.01,4))
  # windows.x <- seq(min(rf$value),max(rf$value),length.out=11)
  windows.x <- quantile(rf$value,probs=seq(0,1,0.025))
  x.avg <- (windows.x[-1] + windows.x[-length(windows.x)])/2
  cols <- viridis_pal()(5)
  r <- range(windows.x)*c(1/1.01,1.01)
  
  x <- sort(rep(windows.x,2))
  x <- x[-length(x)]
  x <- x[-1]
  
  # GY93
  cat("GY93 max power:\n")
  plot(NULL,NULL,xlim=r,ylim=c(0,1),xlab="posterior mean RF",ylab="power",main="GY93")
  for (i in 1:5) {
    d <- all.d[i]
    
    pvals <- read.csv("~/Downloads/2020-01-06_csv/agg_G93_pvalue.csv",stringsAsFactors=FALSE)
    pvals <- pvals[pvals$X.d. == d,]
    detect <- pvals$value < 0.05
    rf <- read.csv("~/Downloads/2020-01-06_csv/agg_rf.mean.csv",stringsAsFactors=FALSE)
    rf <- rf$value[rf$X.d. == d]
    
    windows.y <- sapply(2:length(windows.x),function(i){
      xl <- windows.x[i-1]
      xh <- ifelse(i == length(windows.x),1000,windows.x[i])
      indices <- rf < xh & rf >= xl
      sum(detect[indices])/sum(indices)
    })
    cat("  ",max(windows.y),"\n")
    # lines(x.avg,windows.y,col=cols[i],lwd=2)
    
    y <- as.numeric(rbind(windows.y,windows.y))
    
    lines(x,y,col=cols[i],lwd=2)
    
  }
  legend("topleft",legend=rev(paste0("d = ",all.d)),fill=rev(cols),bty="n",border=NA,cex=1)
  
  # kurt
  cat("kurt max power:\n")
  plot(NULL,NULL,xlim=r,ylim=c(0,1),xlab="posterior mean RF",ylab="power",main="kurtosis(MI)")
  for (i in 1:5) {
    d <- all.d[i]
    
    pvals <- read.csv("~/Downloads/2020-01-06_csv/agg_kurtosis_pvalue.csv",stringsAsFactors=FALSE)
    pvals <- pvals[pvals$X.d. == d,]
    detect <- pvals$value < 0.05
    rf <- read.csv("~/Downloads/2020-01-06_csv/agg_rf.mean.csv",stringsAsFactors=FALSE)
    rf <- rf$value[rf$X.d. == d]
    
    windows.y <- sapply(2:length(windows.x),function(i){
      xl <- windows.x[i-1]
      xh <- ifelse(i == length(windows.x),1000,windows.x[i])
      indices <- rf < xh & rf >= xl
      sum(detect[indices])/sum(indices)
    })
    cat("  ",max(windows.y),"\n")
    
    # lines(x.avg,windows.y,col=cols[i],lwd=2)
    
    y <- as.numeric(rbind(windows.y,windows.y))
    
    lines(x,y,col=cols[i],lwd=2)
    
  }
  
  # max
  cat("max max power:\n")
  plot(NULL,NULL,xlim=r,ylim=c(0,1),xlab="posterior mean RF",ylab="power",main="max(MI)")
  for (i in 1:5) {
    d <- all.d[i]
    
    pvals <- read.csv("~/Downloads/2020-01-06_csv/agg_max_pvalue.csv",stringsAsFactors=FALSE)
    
    pvals <- pvals[pvals$X.d. == d,]
    detect <- pvals$value < 0.05
    rf <- read.csv("~/Downloads/2020-01-06_csv/agg_rf.mean.csv",stringsAsFactors=FALSE)
    rf <- rf$value[rf$X.d. == d]
    
    windows.y <- sapply(2:length(windows.x),function(i){
      xl <- windows.x[i-1]
      xh <- ifelse(i == length(windows.x),1000,windows.x[i])
      indices <- rf < xh & rf >= xl
      sum(detect[indices])/sum(indices)
    })
    cat("  ",max(windows.y),"\n")
    
    # lines(x.avg,windows.y,col=cols[i],lwd=2)
    
    y <- as.numeric(rbind(windows.y,windows.y))
    
    lines(x,y,col=cols[i],lwd=2)
    
  }
dev.off()

