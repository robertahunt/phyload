## First we write all the Rev scripts we'll need

# Script to pool alignments
source("simulation_study/preliminary/0_setup/src/simulation/merge_alignments.R")

# Get Rev source script
revscript <- scan("simulation_study/preliminary/0_setup/src/simulation/rev_model_template.Rev",sep="\n",what="character",strip.white=FALSE)

# Make it a length 1 character vector
revscript <- paste0(revscript,collapse="\n")

# Things that will be the same for all simulations
revscript <- gsub("<<NSITES>>","1496",revscript)

this.seed <- 1233

## Simulate alignments of decreasing site counts, no epistasis
prop.epi <- c(0,0.25,0.5,0.75)
this.d <- 0

for (i in 1:4) {
  this.prop.epi <- prop.epi[i]
  this.seed <- this.seed + 1
  # Make sure we have somewhere to write to
  this.dir <- paste0("simulation_study/preliminary/1_decreasing_site_counts/frac-",this.prop.epi)
  if (!file.exists(this.dir)) {
    dir.create(this.dir)
    dir.create(paste0(this.dir,"/src"))
    dir.create(paste0(this.dir,"/data"))
    dir.create(paste0(this.dir,"/output"))
  }
  # Fill out the rev template for this simulation cell
  this.revscript <- revscript
  this.revscript <- gsub("<<SEED>>",this.seed,this.revscript)
  this.revscript <- gsub("<<PROP_EPI>>",this.prop.epi,this.revscript)
  this.revscript <- gsub("<<EPISTASIS_D>>",this.d,this.revscript)
  this.revscript <- gsub("<<TARGET_DIRECTORY>>",this.dir,this.revscript)
  
  cat(this.revscript,file=paste0(this.dir,"/src/simulate.Rev"))
  
  # Script simulates two components of the alignment, the site-IID sites and epistatic sites
  system2(command="rb",args=paste0(this.dir,"/src/simulate.Rev"))

  # Put the alignments into the one we'll analyze
  for (rep in 1:100) {
    system2("mv",args=c(paste0(this.dir,"/data/base_alignment_",rep,".nex"),paste0(this.dir,"/data/",rep,".nex")))
  }

  # Remove the unneeded alignments
  system2("rm",paste0(this.dir,"/data/epistatic_alignment_*"))
  cat(this.dir,"\n")
  
}


## Simulate the grid of d and proportion of epistatic sites
this.seed <- 8471
prop.epi <- c(0.25,0.5,0.75,1)
d.vals <- c(0,1,2,5,10)

for (i in 1:4) {
  this.prop.epi <- prop.epi[i]
  for (j in 1:5) {
    this.d <- d.vals[j]
    this.seed <- this.seed + 1
    # Make sure we have somewhere to write to
    this.dir <- paste0("simulation_study/preliminary/2_epistasis/prop-",this.prop.epi,"_d-",this.d)
    if (!file.exists(this.dir)) {
      dir.create(this.dir)
      dir.create(paste0(this.dir,"/src"))
      dir.create(paste0(this.dir,"/data"))
      dir.create(paste0(this.dir,"/output"))
    }
    # Fill out the rev template for this simulation cell
    this.revscript <- revscript
    this.revscript <- gsub("<<SEED>>",this.seed,this.revscript)
    this.revscript <- gsub("<<PROP_EPI>>",this.prop.epi,this.revscript)
    this.revscript <- gsub("<<EPISTASIS_D>>",this.d,this.revscript)
    this.revscript <- gsub("<<TARGET_DIRECTORY>>",this.dir,this.revscript)

    cat(this.revscript,file=paste0(this.dir,"/src/simulate.Rev"))

    # Simulate the two components of the alignment, the site-IID sites and epistatic sites
    system2(command="rb",args=paste0(this.dir,"/src/simulate.Rev"))

    # Put the alignments into the one we'll analyze
    for (rep in 1:100) {
      this.base <- paste0(this.dir,"/data/base_alignment_",rep,".nex")
      this.epi <- paste0(this.dir,"/data/epistatic_alignment_",rep,".nex")
      mergeAlignments(this.base,this.epi,paste0(this.dir,"/data/",rep,".nex"))
    }

    # Remove the unneeded alignments
    system2("rm",paste0(this.dir,"/data/base_alignment_*"))
    system2("rm",paste0(this.dir,"/data/epistatic_alignment_*"))

    cat(this.dir,"\n")

  }
}
