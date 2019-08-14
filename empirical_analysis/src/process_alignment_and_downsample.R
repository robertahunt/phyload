library(phangorn)

set.seed(47)

# Function to convert AA,AC,...,TT to 0,1,...,A,B,...,F
# If either first or second base are ambiguous, returns "-"
convertToDoublet <- function(i,j) {
  i <- tolower(i)
  j <- tolower(j)
  if (i == "a") {
    if (j == "a") {
      return("0")
    } else if (j == "c") {
      return("1")
    } else if (j == "g") {
      return("2")
    } else if (j == "t") {
      return("3")
    } else {
      return("-")
    }
  } else if (i == "c") {
    if (j == "a") {
      return("4")
    } else if (j == "c") {
      return("5")
    } else if (j == "g") {
      return("6")
    } else if (j == "t") {
      return("7")
    } else {
      return("-")
    }
  } else if (i == "g") {
    if (j == "a") {
      return("8")
    } else if (j == "c") {
      return("9")
    } else if (j == "g") {
      return("A")
    } else if (j == "t") {
      return("B")
    } else {
      return("-")
    }
  } else if (i == "t") {
    if (j == "a") {
      return("C")
    } else if (j == "c") {
      return("D")
    } else if (j == "g") {
      return("E")
    } else if (j == "t") {
      return("F")
    } else {
      return("-")
    }
  } else {
    return("-")
  }
}

doublet.chars <- c(as.character(0:9),"A","B","C","D","E","F")

# Get original alignment as nexus alignment
aln <- read.phyDat("empirical_analysis/data/Tsagkogeorga-BMCEvolBiol2010-Tunicata_18S_110taxa.nex","nexus")
original.nsites <- dim(as.character(aln))[2]

# Get original alignment file as character vector (to get character sets)
aln.strings <- scan("empirical_analysis/data/Tsagkogeorga-BMCEvolBiol2010-Tunicata_18S_110taxa.nex",what=character(),sep="\n")

all.ends <- grep(";",aln.strings)

# Information on sites we need
exclude.first <- grep("exclude",aln.strings)
exclude.last <- all.ends[min(which(all.ends > exclude.first))]
exclude <- aln.strings[exclude.first:exclude.last]
exclude <- gsub("[exclude ","",exclude,fixed=TRUE)
exclude <- gsub(";]","",exclude,fixed=TRUE)
exclude <- unlist(lapply(exclude,function(line){
  vec <- strsplit(line," ")[[1]]
  unlist(lapply(vec,function(x){
    if (grepl("-",x)) {
      x <- as.numeric(strsplit(x,"-",fixed=TRUE)[[1]])
      return(x[1]:x[2])
    } else {
      return(as.numeric(x))
    }
  }))
}))

loop.first <- grep("LOOP",aln.strings)
loop.last <- all.ends[min(which(all.ends > loop.first))]
loop <- aln.strings[loop.first:loop.last]
loop <- gsub("[CHARSET  LOOP = ","",loop,fixed=TRUE)
loop <- gsub(";]","",loop,fixed=TRUE)
loop <- unlist(lapply(loop,function(line){
  vec <- strsplit(line," ")[[1]]
  unlist(lapply(vec,function(x){
    if (grepl("-",x)) {
      x <- as.numeric(strsplit(x,"-",fixed=TRUE)[[1]])
      return(x[1]:x[2])
    } else {
      return(as.numeric(x))
    }
  }))
}))

pair.first <- grep("PAIRES",aln.strings)
pair.last <- all.ends[min(which(all.ends > pair.first))]
pair <- aln.strings[pair.first:pair.last]
pair <- gsub("[PAIRES SITES ","",pair,fixed=TRUE)
pair <- gsub(";]","",pair,fixed=TRUE)
pair <- lapply(pair,function(line){
  vec <- strsplit(line,", ")[[1]]
  vec <- gsub(",","",vec)
  sapply(vec,function(x){
    as.numeric(strsplit(x,":",fixed=TRUE)[[1]])
  })
})
pair <- do.call(cbind,pair)

# Remove outgroups
outgroups <- c("HomoXsapie","XenopusXla","DanioXreri","Strongyloc","Saccogloss",
               "Ptychodera","Branchiost","Petromyzon","ChrysemysX","RajaXschmi",
               "GallusXgal","AnolisXcar","AntedonXse","Balanoglos","AsteriasXa")

aln <- as.character(aln)

aln <- aln[!row.names(aln) %in% outgroups,]

# Downsample to 50 taxa
aln <- aln[sample.int(95,50),]

# Extract loop region, write alignment
aln.loop <- aln[,c((c(1:original.nsites) %in% loop) & !c(1:original.nsites) %in% exclude)]
write.phyDat(as.phyDat(aln.loop),"empirical_analysis/data/loop_50.nex","nexus")

# Extract pairs and recode, write alignment
aln.pair <- matrix(nrow=dim(aln)[1],ncol=dim(pair)[2])
row.names(aln.pair) <- row.names(aln)
for (i in 1:dim(aln.pair)[1]) {
  for (j in 1:dim(aln.pair)[2]) {
    x <- aln[i,pair[1,j]]
    y <- aln[i,pair[2,j]]
    aln.pair[i,j] <- convertToDoublet(x,y)
  }
}
write.phyDat(phyDat(aln.pair,type="USER",levels=doublet.chars,ambiguity="-"),"empirical_analysis/data/pair_50.nex","nexus")

# Correct the datatype
aln.pair <- scan("empirical_analysis/data/pair_50.nex",what=character(),sep="\n",strip.white=FALSE)
aln.pair <- gsub("PROTEIN","Standard",aln.pair)
aln.pair <- gsub("INTERLEAVE=YES","INTERLEAVE=YES SYMBOLS=\"0123456789ABCDEF\"",aln.pair)
cat(aln.pair,sep="\n",file="empirical_analysis/data/pair_50.nex")

# Write reduced all-nucleotide alignment as fasta
write.phyDat(phyDat(aln),"empirical_analysis/data/unpartitioned_50.fasta",format="fasta")

