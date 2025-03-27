#!/usr/bin/env Rscript

###############################################################################
# 
###############################################################################
###############################################################################
# Set parameters

QUAL_thresh <- 1000
triploid_Lines <- c("F_B", "F_E")

no_Anc <- c("H_F", "H_H", "N_F", "N_H")

contaminated <- c("N_B08", "N_B12", "N_D08", "N_D10", "N_E09", "N_E10", "N_G03", 
                  "H_D05", "H_B12", "H_E03", "H_G05", "H_G07")

bad_seq <- c("N_A01", "N_B12", "N_C03")

# End Section #################################################################
###############################################################################
#                                                                             #
#                                                                             #
###############################################################################
###############################################################################
# Load in required functions and packages

lib_dir <- "~/R/x86_64-pc-linux-gnu-library/3.6"
## Load packages
library(optparse, lib.loc = lib_dir)
library(stringr, lib.loc = lib_dir)
library(rlang, lib.loc = lib_dir)
library(vctrs, lib.loc = lib_dir)
library(cli, lib.loc = lib_dir)
library(dplyr, lib.loc = lib_dir)

# End Section #################################################################
###############################################################################
#                                                                             #
#                                                                             #
###############################################################################
###############################################################################
# Get command line arguements

# Create set of input arguements to parse
option_list = list(
    make_option(c("-s", "--SNPsIn"), type = "character", default = NULL,
                action = "store", help = "SNPs table input"),
    make_option(c("-S", "--SNPsOut"), type = "character", default = NULL, 
                action = "store", help = "output SNP table file (must be .RData or .csv)")
) 

# Parse arguements and create list of values
opt_parser <- OptionParser(option_list = option_list, add_help_option = FALSE)
input_args <- parse_args(opt_parser)

# Determine whether arguements are valid and stop script if not
args_invalid <- any(sapply(input_args, function(x) is.null(x)))
if(args_invalid) {
    print("One or more options not provided")
} else {

# Assign aguement list entries to separate variables
# SNPs_input <- "/oasis/tscc/scratch/mioverto/LOH_methods/Overton_etal/variants/final/all_BYxRM_reconciled_SNPs.RData"
# SNPs_output <- "/oasis/tscc/scratch/mioverto/LOH_methods/Overton_etal/variants/final/all_BYxRM_clean_SNPs.RData"
SNPs_input <- as.character(input_args[1])
SNPs_output <- as.character(input_args[2])

# End Section #################################################################
###############################################################################
#                                                                             #
#                                                                             #
###############################################################################
###############################################################################

input_RData <- grepl(".RData", SNPs_input)
if(input_RData) {
    load(file = SNPs_input)
} else {
   SNPs_merge_finalGT <- read.table(file = SNPs_input, header = T, sep = ",")
}

SNPs_merge_finalGT <- SNPs_merge_finalGT %>% 
  ungroup() %>%
  mutate(Line = substr(ID, 1, 3)) %>%
  filter(!Line %in% c(triploid_Lines), !ID %in% c(contaminated, bad_seq)) %>%
  filter(QUAL_P1call >= QUAL_thresh | QUAL_P2call >= QUAL_thresh)

output_RData <- grepl(".RData", SNPs_output)
if(output_RData) {
    save(SNPs_merge_finalGT, file = SNPs_output)
} else {
   write.table(SNPs_merge_finalGT, file = SNPs_output, sep = ",",
             col.names = T, row.names = F)
}

}