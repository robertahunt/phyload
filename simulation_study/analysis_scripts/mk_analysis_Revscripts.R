# This is designed to be called with Rscript with several arguments
# arg1: number of scripts to make
# arg2: seed
# arg3: path to output alignments (relative to location of phyload repository)

# Call this script from phyload/simulation_study

# Get arguments
args = commandArgs(trailingOnly=TRUE)

if (!length(args) == 5) {
  stop("This script requires 5 arguments")
}

nscripts  <- args[1]
this.seed <- args[2]
out.file  <- args[3]

# Directory that we're outputting our Revscript to (to pass to Revscript)
out.dir <- dirname(out.file)

# Get Rev source script
revscript <- scan("simulation_scripts/analysis_template.Rev",sep="\n",what="character",strip.white=FALSE)

# Make it a length 1 character vector
revscript <- paste0(revscript,collapse="\n")

for (i in 1:nscripts) {
  
  # Fill out the rev template for this simulation cell
  this.revscript <- revscript
  this.revscript <- gsub("<<SEED>>",this.seed,this.revscript)
  this.revscript <- gsub("<<REPLICATE>>",i,this.revscript)
  
  cat(this.revscript,file=paste0(out.dir,"/analyze_",i,".Rev"))
  
  this.seed <- this.seed + 1
}
