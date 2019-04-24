# This is designed to be called with Rscript with several arguments
# arg1: number of total sites
# arg2: proportion of epistatic sites
# arg3: d parameter for model
# arg4: path to output alignments (relative to location of phyload repository)
# arg5: path to installed version of RevBayes (can be "rb" if installed)
# arg6: seed

# Call this script from the top level of the phyload repository

# Get arguments
args = commandArgs(trailingOnly=TRUE)

if (!length(args) == 6) {
  stop("This script requires 6 arguments")
}

nsites    <- args[1]
prop.epi  <- args[2]
this.d    <- args[3]
out.dir   <- args[4]
rb.path   <- args[5]
this.seed <- args[6]


# Script to pool alignments
source("simulation_study/simulation_scripts/merge_alignments.R")

# Get Rev source script
revscript <- scan("simulation_study/simulation_scripts/rev_model_template.Rev",sep="\n",what="character",strip.white=FALSE)

# Make it a length 1 character vector
revscript <- paste0(revscript,collapse="\n")

# Things that will be the same for all simulations
revscript <- gsub("<<NSITES>>",as.character(nsites),revscript)

# Make sure we have somewhere to write to
if (!file.exists(out.dir)) {
  dir.create(out.dir)
}

# Fill out the rev template for this simulation cell
this.revscript <- revscript
this.revscript <- gsub("<<SEED>>",this.seed,this.revscript)
this.revscript <- gsub("<<PROP_EPI>>",prop.epi,this.revscript)
this.revscript <- gsub("<<EPISTASIS_D>>",this.d,this.revscript)
this.revscript <- gsub("<<TARGET_DIRECTORY>>",out.dir,this.revscript)

cat(this.revscript,file=paste0(out.dir,"/simulate_alignments.Rev"))

# Script simulates two components of the alignment, the site-IID sites and epistatic sites
system2(command=as.character(rb.path),args=paste0(out.dir,"/simulate_alignments.Rev"))

# Put the alignments into the one we'll analyze
for (rep in 1:100) {
  system2("mv",args=c(paste0(out.dir,"/base_alignment_",rep,".nex"),paste0(out.dir,"/",rep,".nex")))
}

# Remove the unneeded alignments
system2("rm",paste0(out.dir,"/epistatic_alignment_*"))
