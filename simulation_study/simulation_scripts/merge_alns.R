# This is designed to be called with Rscript with several arguments
# arg1: alignment 1 path
# arg2: alignment 2 path
# arg3: output alignment path

# This script concatenates two nexus alignments and outputs the concatenated
# alignment

# Get arguments
args = commandArgs(trailingOnly=TRUE)

if (!length(args) == 3) {
  stop("This script requires 2 arguments")
}

# Function for merging
mergeAlignments <- function(regular, epistatic, output) {
  # Get both alignments
  reg <- scan(regular,what=character(),sep="\n")
  epi <- scan(epistatic,what=character(),sep="\n")

  # recover()

  # Check for empty alignments
  reg_empty <- any(grepl("nchar=0",tolower(reg)))
  epi_empty <- any(grepl("nchar=0",tolower(epi)))
  
  if ( !epi_empty ) {
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
      epi_nchar <- 2 * length(doublets)
      epi_data <- epi[epi_aln_lines]
    }
  } else {
    epi_nchar <- 0
    epi_data <- ""
  }

  # Put them together (shove epistasis alignment lines after regular ones), amend number of characters
  reg_aln_start <- grep("Matrix",reg)+1
  reg_aln_end <- grep(";",reg)
  reg_aln_end <- min(reg_aln_end[reg_aln_end > reg_aln_start]) - 1
  reg_last_lines <- reg[(reg_aln_end + 1):length(reg)]
  reg_data <- reg[reg_aln_start:reg_aln_end]

  merged <- reg[1:(reg_aln_start-1)]
  
  if ( reg_empty && epi_empty ) {
    # Empty alignmnent should still have taxon names
    merged <- c(merged,reg_data)
  } else {
    if ( !reg_empty ) {
      merged <- c(merged,reg_data)
    }
    if ( !epi_empty ) {
      merged <- c(merged,epi_data)
    }
  }
  
  merged <- c(merged,reg_last_lines)
  
  merged_char_line <- grep("nchar",merged)
  merged_nchar <- merged[merged_char_line]
  merged_nchar <- strsplit(merged_nchar,"nchar=")[[1]]
  total_nchar <- as.numeric(gsub(";","",merged_nchar[2])) + epi_nchar
  total_nchar <- paste0(merged_nchar[1],"nchar=",total_nchar,";")
  merged[merged_char_line] <- total_nchar

  # Length 0 alignments produce no datatype, we must add it if it went missing
  merged_format_line <- grep("format",tolower(merged))
  if ( grepl("datatype= ",merged[merged_format_line]) ) {
    merged[merged_format_line] <- gsub("datatype= ","datatype=DNA ",merged[merged_format_line])
  }

  merged <- gsub("Format ","Format interleave=yes ",merged)

  cat(merged, file=output, sep="\n")

}

# Call merge function
mergeAlignments(args[1], args[2], args[3])
