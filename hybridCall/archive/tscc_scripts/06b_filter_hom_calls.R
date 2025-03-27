#!/usr/bin/env Rscript

###############################################################################
# This script uses the parameters generated from 06a_train_het_filter.R to 
# calculate the Alt allele depth thresholds for total depth to reject
# heterozygous calls from the SNP table output of either 
# 05_reconcile_callsets.R or remove_prob_clones.R
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

# This script uses frunctions from the 05a_reconcile_functions.R script, which
# must be located in the same folder

lib_dir <- "~/R/x86_64-pc-linux-gnu-library/3.6"
## Load packages
library(optparse, lib.loc = lib_dir)
library(stringr, lib.loc = lib_dir)
library(rlang, lib.loc = lib_dir)
library(emdbook, lib.loc = lib_dir)
library(vctrs, lib.loc = lib_dir)
library(cli, lib.loc = lib_dir)
library(dplyr, lib.loc = lib_dir, warn.conflicts = F)


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
# Get command line arguements #################################################

# Create set of input arguements to parse
option_list <- list(
    make_option(c("-s", "--SNPsInput"), type = "character", default = NULL,
                help = "SNPs table input", metavar = "character"),
    make_option(c("-m", "--modelInput"), type = "character", default = NULL,
                help = "Het threshold model table", metavar = "character"),
    make_option(c("-r", "--FPrate"), type = "character", default = NULL,
                help = "False positive rate under model (Rate of exclusion of true Hets)", metavar = "numeric"),
    make_option(c("-S", "--SNPsOutput"), type = "character", default = NULL,
                help = "SNPs table output", metavar = "character")
    
) 

# Parse arguements and create list of values
opt_parser <- OptionParser(option_list = option_list)
input_args <- parse_args(opt_parser)

# Determine whether arguements are valid and stop script if not
args_invalid <- any(sapply(input_args, function(x) is.null(x)))
if(args_invalid) {
    print("One or more options not provided")
} else {

# Assign aguement list entries to separate variables
# SNPs_input <- "/oasis/tscc/scratch/mioverto/LOH_methods/Overton_etal/variants/reconciled/all_BYxRM_hetFilter_SNPs.RData"
# model_input <- "/oasis/tscc/scratch/mioverto/LOH_methods/Overton_etal/metadata/hom_model.RData"
# FE_rate <- 0.05
# SNPs_output <- "/oasis/tscc/scratch/mioverto/LOH_methods/Overton_etal/variants/reconciled/all_BYxRM_hetFilter_SNPs.RData"
SNPs_input <- as.character(input_args[1])
model_input <- as.character(input_args[2])
FE_rate <- as.numeric(input_args[3])
SNPs_output <- as.character(input_args[4])

# End Section #################################################################
###############################################################################
#                                                                             #
#                                                                             #
###############################################################################
# Load in SNP data frame and model list #######################################

input_RData <- grepl(".RData", SNPs_input)
if(input_RData) {
    load(file = SNPs_input)
} else {
   SNPs_merge_finalGT <- read.table(file = SNPs_input, header = T, sep = ",")
}

# model input should be a .RData list with two elements:
# 1 - a vector of model parameters: lin_p and exp_p
# 2 - a dataframe of depth bins and inferred Lambda values
load(file = model_input)
hom_params <- hom_model_output[[1]]
model_table <- hom_model_output[[2]]

# End Section #################################################################
###############################################################################
#                                                                             #
#                                                                             #
###############################################################################
# Calculate model thresholds ###########################

max_dp <- SNPs_merge_finalGT %>% filter(GT %in% c("0/0", "1/1")) %>% pull(Sum_DP_final) %>% max()

Alt_DP_thrsh <- data.frame(Sum_DP = 1:max_dp) %>% 
  mutate(thrsh = LambdaLMthrsh(Sum_DP, hom_params[1], hom_params[2], .alpha = (1 - FE_rate)))

# End Section #################################################################
###############################################################################
#                                                                             #
#                                                                             #
###############################################################################
# Load in SNP data frame and apply model thresholds ###########################

SNPs_merge_finalGT$Pois_reject <- F
for(th in unique(Alt_DP_thrsh$thrsh)) {
    # th <- 1
    min_dp <- Alt_DP_thrsh %>% filter(thrsh == th) %>% pull(Sum_DP) %>% min(.)
    max_dp <- Alt_DP_thrsh %>% filter(thrsh == th) %>% pull(Sum_DP) %>% max(.)
    i_dp_R <- SNPs_merge_finalGT$GT == "0/0" &
            SNPs_merge_finalGT$Sum_DP_final >= min_dp & 
            SNPs_merge_finalGT$Sum_DP_final <= max_dp & 
            SNPs_merge_finalGT$Alt_DP_final > th
    i_dp_A <- SNPs_merge_finalGT$GT == "1/1" &
            SNPs_merge_finalGT$Sum_DP_final >= min_dp & 
            SNPs_merge_finalGT$Sum_DP_final <= max_dp & 
            SNPs_merge_finalGT$Ref_DP_final > th
    SNPs_merge_finalGT$Pois_reject[i_dp_R | i_dp_A] <- T
  
}

output_RData <- grepl(".RData", SNPs_output)
if(output_RData) {
    save(SNPs_merge_finalGT, file = SNPs_output)
} else {
   write.table(SNPs_merge_finalGT, file = SNPs_output, sep = ",",
             col.names = T, row.names = F)
}

}