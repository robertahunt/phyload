# This is designed to be called with Rscript with several arguments
# arg1: the (relative) file path to directory we're simulating in
# arg2: do we want to delete the component alignments?

# Call this script from phyload/simulation_study

# Get arguments
args = commandArgs(trailingOnly=TRUE)

if (!length(args) >= 1) {
  stop("This script requires at least 1 argument")
}

out.dir <- args[1]

if (length(args) == 2) {
  delete <- args[2]  
} else {
  delete <- FALSE
}


replicate <- basename(dirname(out.dir))

iid.aln  <- paste0(out.dir,"/iid_aln.nex")
epi.aln  <- paste0(out.dir,"/epi_aln.nex")
out.file <- paste0(out.dir,"/aln.nex")

# Function for merging
mergeAlignments <- function(regular,epistatic,out.file) {
  # Get both alignments
  reg <- scan(regular,what=character(),sep="\n")
  epi <- scan(epistatic,what=character(),sep="\n")
  
  # recover()
  
  # Turn Rev's "standard" datatype into doublets
  epi_aln_start <- grep("Matrix",epi)+1
  epi_aln_end <- grep(";",epi)
  epi_aln_end <- min(epi_aln_end[epi_aln_end > epi_aln_start]) - 1
  epi_aln_lines <- epi_aln_start:epi_aln_end
  
  for (i in epi_aln_lines) {
    tmp <- strsplit(epi[i]," ")[[1]]
    nspaces <- length(grep("^$",tmp)) + 1 # We split on space, so between two spaces is a blank
    taxon <- tmp[1]
    doublets <- strsplit(tmp[4],"")[[1]]
    for (j in 1:length(doublets)) {
      if ( doublets[j] == "0" ) {
        doublets[j] <- "AA"
      } else if ( doublets[j] == "1" ) {
        doublets[j] <- "AC"
      } else if ( doublets[j] == "2" ) {
        doublets[j] <- "AG"
      } else if ( doublets[j] == "3" ) {
        doublets[j] <- "AT"
      } else if ( doublets[j] == "4" ) {
        doublets[j] <- "CA"
      } else if ( doublets[j] == "5" ) {
        doublets[j] <- "CC"
      } else if ( doublets[j] == "6" ) {
        doublets[j] <- "CG"
      } else if ( doublets[j] == "7" ) {
        doublets[j] <- "CT"
      } else if ( doublets[j] == "8" ) {
        doublets[j] <- "GA"
      } else if ( doublets[j] == "9" ) {
        doublets[j] <- "GC"
      } else if ( doublets[j] == "A" ) {
        doublets[j] <- "GG"
      } else if ( doublets[j] == "B" ) {
        doublets[j] <- "GT"
      } else if ( doublets[j] == "C" ) {
        doublets[j] <- "TA"
      } else if ( doublets[j] == "D" ) {
        doublets[j] <- "TC"
      } else if ( doublets[j] == "E" ) {
        doublets[j] <- "TG"
      } else if ( doublets[j] == "F" ) {
        doublets[j] <- "TT"
      } 
    }
    epi[i] <- paste0(taxon,paste0(rep(" ",nspaces),collapse=""),paste0(doublets,collapse=""))
  }
  
  # Put them together (shove epistasis alignment lines after regular ones), amend number of characters
  # Turn Rev's "standard" datatype into doublets
  reg_aln_start <- grep("Matrix",reg)+1
  reg_aln_end <- grep(";",reg)
  reg_aln_end <- min(reg_aln_end[reg_aln_end > reg_aln_start]) - 1
  reg_last_lines <- reg[(reg_aln_end + 1):length(reg)]
  reg <- reg[1:reg_aln_end]
  reg <- c(reg,epi[epi_aln_lines],reg_last_lines)
  
  reg_char_line <- grep("nchar",reg)
  reg_nchar <- reg[reg_char_line]
  reg_nchar <- strsplit(reg_nchar,"nchar=")[[1]]
  total_nchar <- as.numeric(gsub(";","",reg_nchar[2])) + 2 * length(doublets)
  total_nchar <- paste0(reg_nchar[1],"nchar=",total_nchar,";")
  reg[reg_char_line] <- total_nchar
  
  reg <- gsub("Format ","Format interleave=yes ",reg)
  
  cat(reg,file=out.file,sep="\n")
  
}

# Call merge function
mergeAlignments(iid.aln,epi.aln,out.file)

if (delete) {
  system2("rm",iid.aln)
  system2("rm",epi.aln)
}