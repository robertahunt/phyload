# This is designed to be called with Rscript with several arguments
# arg1: number of scripts to make
# arg2: seed
# arg3: path to directory to output the revscript
# arg4: path to rb to call it

# Call this script from phyload/simulation_study

# Get arguments
args = commandArgs(trailingOnly=TRUE)

if (!length(args) == 4) {
  stop("This script requires 4 arguments")
}

nscripts  <- args[1]
this.seed <- args[2]
out.dir   <- args[3]
rb.path   <- args[4]


# Get Rev source script
revscript <- scan("simulation_scripts/analysis_template.Rev",sep="\n",what="character",strip.white=FALSE)

# Make it a length 1 character vector
revscript <- paste0(revscript,collapse="\n")

# Print revscripts
for (i in 1:nscripts) {
  
  # Fill out the rev template for this simulation cell
  this.revscript <- revscript
  this.revscript <- gsub("<<SEED>>",this.seed,this.revscript)
  this.revscript <- gsub("<<REPLICATE>>",i,this.revscript)
  
  cat(this.revscript,file=paste0(out.dir,"/analyze_",i,".Rev"))
  
  this.seed <- this.seed + 1
}

# Run revscripts
for (i in 1:nscripts) {
  system2(command=as.character(rb.path),args=paste0(out.dir,"/analyze_",i,".Rev"))
}

