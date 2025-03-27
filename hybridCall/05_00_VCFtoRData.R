#!/usr/bin/env Rscript

###############################################################################
# This script loads in a VCF file and parses it into a data.frame.
# The expected VCF format is a multi-sample VCF output from GATK, and may 
# contain a combination of founders and/or end-point clones.
# The output is a single data.frame containing all of the rows in the VCF.
###############################################################################
###############################################################################
# Set parameters


# End Section #################################################################
###############################################################################
#                                                                             #
#                                                                             #
###############################################################################
###############################################################################
# Load in required functions and packages

lib_dir <- "/tscc/nfs/home/mioverto/bin/R/x86_64-pc-linux-gnu-library/4.1"
.libPaths(lib_dir)

## Load packages
library(optparse)
library(stringr)
library(rlang)
library(vctrs)
library(cli)
library(tidyr)
library(dplyr)


# This script uses frunctions from the 05a_reconcile_functions.R script, which
# must be located in the same folder
code_dir <- "~/code/hybridCall/"
source(paste0(code_dir, "05a_reconcile_functions.R"))

# End Section #################################################################
###############################################################################
#                                                                             #
#                                                                             #
###############################################################################
###############################################################################
# Get command line arguements

# Create set of input arguements to parse
option_list = list(
    make_option(c("-i", "--VCFin"), type = "character", default = NULL,
                action = "store", help = "VCF file name"),
    make_option(c("-o", "--RDataOut"), type = "character", default = NULL,
                action = "store", help = "RData output file name")
) 

# Parse arguements and create list of values
opt_parser <- OptionParser(option_list = option_list, add_help_option = FALSE)
input_args <- parse_args(opt_parser)

# Assign aguement list entries to separate variables
# proj_dir <- "/tscc/lustre/ddn/scratch/mioverto/LOH_methods/Overton_etal/"
# VCF_in <- paste0(proj_dir, "variants/final/BY/all_BY.filter.vcf")
# VCF_in <- paste0(proj_dir, "variants/final/RM/all_RM.filter.vcf")
# RData_out <- paste0(proj_dir, "variants/final/BY/all_BY.RData")
# RData_out <- paste0(proj_dir, "variants/final/RM/all_RM.RData")

VCF_in <- as.character(input_args[1])
RData_out <- as.character(input_args[2])

# Determine whether arguements are valid and stop script if not
args_invalid <- any(sapply(input_args, function(x) is.null(x)))
if(args_invalid) {
    print("One or more options not provided")
} else {

# Get parent ID
VCF_filename <- str_split(VCF_in, "") %>% unlist()
i_undscr <- grep("_", VCF_filename)
i_dot <- grep('\\.', VCF_filename)
P1 <- VCF_filename[(i_undscr[length(i_undscr)] + 1):(i_dot[1] - 1)] %>% paste(collapse = "")

# End Section #################################################################
###############################################################################
#                                                                             #
#                                                                             #
###############################################################################
###############################################################################
# Import VCF file

# Load in VCF file
VCF_df <- HC_multiVcfParse_allVar(VCF_in)
VCF_df$ID <- gsub(paste0("_", P1), "", VCF_df$ID)

cat("VCF converted to data.frame\n")

# Save data.frame
save(VCF_df, file = RData_out)

cat("RData file saved\n")