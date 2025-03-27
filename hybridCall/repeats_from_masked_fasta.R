#!/usr/bin/env Rscript

###############################################################################
# This script loads in a soft-masked fasta sequence and extracts the positions 
# of contiguous masked bases.
# The output is a single bed-formatted table of masked positions.
###############################################################################
###############################################################################
# Load and create required functions and libraries ############################

lib_dir <- "/tscc/nfs/home/mioverto/bin/R/x86_64-pc-linux-gnu-library/4.1"
.libPaths(lib_dir)

## Load packages
library(optparse)
library(stringr)
library(rlang)
library(vctrs)
library(cli)
library(dplyr)


## Create functions

# MaskedFastaToBed is the primary function that converts the masked sequence
# to a table of positions
MaskedFastaToBed <- function(fasta_file, format_chrs = F) {
  # fasta_file <- mask_seq
  seqIn <- readLines(fasta_file)
  bed_out <- data.frame(NULL)
  iChrom <- grep(">", seqIn)
  chrom_adj <- 0
  for(i in seq_along(iChrom)) {
    # i <- 1
    chrom_name <- strsplit(seqIn[iChrom[i]], " ")[[1]][1] %>% str_remove(., ">")
    seq_start <- iChrom[i] + 1
    if (i == length(iChrom)) {
      seq_end <- length(seqIn)
    } else {
      seq_end <- (iChrom[i + 1] - 1)
    }
    chrom_seq <- paste0(seqIn[seq_start:seq_end], collapse = "")
    # index all masked bases on chr
    iMasked <- grep("[a-z]", strsplit(chrom_seq, "")[[1]])
    if(length(iMasked) > 0) {
    # within the masked regions, index on the starts of gaps
      iStart_vals <- c(1, (which(diff(iMasked) > 1) + 1))
      # this gives the masked region start indices
      start_vals <- iMasked[iStart_vals] - 1
      start_posi <- start_vals + chrom_adj[i]
      iEnd_vals <- c(which(diff(iMasked) > 1))
      end_vals <- c(iMasked[iEnd_vals], max(iMasked))
      end_posi <- end_vals + chrom_adj[i]
      chrom_df <- data.frame(CHROM = chrom_name, start = start_vals, end = end_vals, 
                            start_POSi = start_posi, end_POSi = end_posi)
      chrom_adj[i + 1] <- nchar(chrom_seq) + chrom_adj[i]
      bed_out <- rbind(bed_out, chrom_df)
    }
  }
  bed_out$CHROM <- factor(bed_out$CHROM)
  bed_out <- bed_out[!is.na(bed_out$start) & !is.infinite(bed_out$end), ]
  if(format_chrs) {
    bed_out$CHROM <- str_pad(as.numeric(bed_out$CHROM), 2, pad = "0")
  }
  return(bed_out)
}

# End Section #################################################################
###############################################################################
#                                                                             #
#                                                                             #
###############################################################################
###############################################################################
# Get command line arguements

# Create set of input arguements to parse
option_list = list(
    make_option(c("-s", "--seq"), type = "character", default = NULL,
                  action = "store", help = "Fasta sequence file"),
    make_option(c("-b", "--bedOut"), type = "character", default = NULL,
                action = "store", help = "Output Bed file"),
    make_option(c("-F", "--formatChrs"), type = "logical",
                action = "store", default = FALSE, help = "Format chromosome names")
) 

# Parse arguements and create list of values
input_args <- parse_args(OptionParser(option_list = option_list))

# Assign aguement list entries to separate variables
mask_seq <- as.character(input_args[1])
output_bed <- as.character(input_args[2])
format_chrs <- as.logical(input_args[3])
# mask_seq <- "/tscc/lustre/ddn/scratch/mioverto/LOH_methods/Pankajam_etal/references/final/S288C/S288C_R64_refseq.fasta"
# output_bed <- "/tscc/lustre/ddn/scratch/mioverto/LOH_methods/Pankajam_etal/references/final/S288C/repeats/S288C_repeats.bed"
# format_chrs <- T

cat(paste0(mask_seq, "\n", output_bed, "\n", ifelse(format_chrs, "True", "False")))

# Determine whether arguements are valid and stop script if not
args_invalid <- any(sapply(input_args, function(x) is.null(x)))
if(args_invalid) {
    print("One or more options not provided")
} else {

# End Section #################################################################
###############################################################################
#                                                                             #
#                                                                             #
###############################################################################
###############################################################################
# Run function on sequence and output BED file

bed_df <- MaskedFastaToBed(mask_seq, format_chrs)

write.table(bed_df, output_bed, sep = "\t", row.names = F, col.names = F, quote = F)

}