# This is designed to be called with Rscript with several arguments
# arg1: seed
# arg2: path to directory to output the revscript
# arg3: path to rb to call it

# Call this script from phyload/simulation_study

# Get arguments
args = commandArgs(trailingOnly=TRUE)

if (!length(args) == 3) {
  stop("This script requires 3 arguments")
}

this.seed <- args[1]
input.aln <- args[2]
out.dir   <- dirname(input.aln)
rb.path   <- args[3]


# Get Rev source script
revscript <- scan("analysis_scripts/analysis_template.Rev",sep="\n",what="character",strip.white=FALSE, quiet=TRUE)

# Make it a length 1 character vector
revscript <- paste0(revscript,collapse="\n")

# Print revscript

# Fill out the rev template for this simulation cell
this.revscript <- revscript
this.revscript <- gsub("<<SEED>>",this.seed,this.revscript)
this.revscript <- gsub("<<TARGET_DIRECTORY>>", out.dir, this.revscript)
this.revscript <- gsub("<<TARGET_ALN>>", input.aln, this.revscript)

cat(this.revscript,file=paste0(out.dir,"/analysis.Rev"))

# Run revscript
system2(command=as.character(rb.path),args=paste0(out.dir,"/analysis.Rev"))
