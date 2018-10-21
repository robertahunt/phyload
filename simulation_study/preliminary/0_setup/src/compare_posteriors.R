library(phangorn)

myblue <- "#4280f4BB"
myorange <- "#ffa035BB"

# The true tree
real <- read.nexus("preliminary_simulating/data/epistasis_simulating_tree.tre")

full.analysis <- read.table("preliminary_simulating/output/full_alignment_analysis.trees",header=TRUE,row.names=1,stringsAsFactors=FALSE)
full.analysis <- lapply(full.analysis$psi,function(treestr){read.tree(text=treestr)})
class(full.analysis) <- "multiPhylo"

regular.analysis <- read.table("preliminary_simulating/output/regular_alignment_analysis.trees",header=TRUE,row.names=1,stringsAsFactors=FALSE)
regular.analysis <- lapply(regular.analysis$psi,function(treestr){read.tree(text=treestr)})
class(regular.analysis) <- "multiPhylo"

full.kf <- KF.dist(full.analysis,real)
regular.kf <- KF.dist(regular.analysis,real)

r <- range(c(full.kf,regular.kf))
breaks <- seq(0.99*r[1],1.01*r[2],length.out=50)

hist(regular.kf,col=myblue,border=NA,breaks=breaks,xlim=r,freq=FALSE,ylim=c(0,150),xlab="KF distance",main="Posterior Distribution of KF distances to true tree")
hist(full.kf,col=myorange,border=NA,breaks=breaks,xlim=r,freq=FALSE,add=TRUE)
legend("topright",legend=c("non-epistatic sites only","full dataset"),fill=c(myblue,myorange),border=NA,bty="n")
