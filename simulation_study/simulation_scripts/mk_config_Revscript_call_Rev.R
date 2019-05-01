# This is designed to be called with Rscript with several arguments
# arg1: number of total sites
# arg2: proportion of epistatic sites
# arg3: d parameter for model
# arg4: seed
# arg5: path to output alignments (relative to location of phyload repository)
# arg6: path to RevBayes

# Call this script from phyload/simulation_study

# Get arguments
args = commandArgs(trailingOnly=TRUE)

if (!length(args) == 6) {
  stop("This script requires 6 arguments")
}

nsites    <- args[1]
prop.epi  <- args[2]
this.d    <- args[3]
this.seed <- args[4]
out.file  <- args[5]
rb.path   <- args[6]

# Directory that we're outputting our Revscript to (to pass to Revscript)
out.dir <- dirname(out.file)

# Get Rev source script
revscript <- scan("simulation_scripts/rev_model_template.Rev",sep="\n",what="character",strip.white=FALSE)

# Make it a length 1 character vector
revscript <- paste0(revscript,collapse="\n")

# Things that will be the same for all simulations
revscript <- gsub("<<NSITES>>",as.character(nsites),revscript)

# Fill out the rev template for this simulation cell
this.revscript <- revscript
this.revscript <- gsub("<<SEED>>",this.seed,this.revscript)
this.revscript <- gsub("<<PROP_EPI>>",prop.epi,this.revscript)
this.revscript <- gsub("<<EPISTASIS_D>>",this.d,this.revscript)
this.revscript <- gsub("<<TARGET_DIRECTORY>>",out.dir,this.revscript)

cat(this.revscript,file=paste0(out.dir,"/simulate_alignments.Rev"))



# Script simulates two components of the alignment, the site-IID sites and epistatic sites
system2(command=as.character(rb.path),args=paste0(out.dir,"/simulate_alignments.Rev"))

