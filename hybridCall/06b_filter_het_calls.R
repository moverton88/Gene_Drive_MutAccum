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

lib_dir <- "/tscc/nfs/home/mioverto/bin/R/x86_64-pc-linux-gnu-library/4.1"
.libPaths(lib_dir)
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
                help = "SNPs table input"),
    make_option(c("-m", "--modelInput"), type = "character", default = NULL,
                help = "Het threshold model table"),
    make_option(c("-t", "--modelType"), type = "character", default = "v",
                help = "Which type of model to use: 
                        'd' for fixed median alpha beta-binomial
                        'n' for fixed mean alpha beta-binomial,
                        or 'v' for beta-binomial from predicted variance under 
                        binomial + error [default]"),
    make_option(c("-r", "--FPrate"), type = "numeric", default = NULL,
                help = "False positive rate under model (Rate of exclusion of 
                        true Hets)"),
    make_option(c("-S", "--SNPsOutput"), type = "character", default = NULL,
                help = "SNPs table output")
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
# SNPs_input <- "/tscc/lustre/ddn/scratch/mioverto/LOH_methods/Loeillet_etal/variants/final/Dual/all_BYxSK1_reconciled_SNPs.RData"
# model_input <- "/tscc/lustre/ddn/scratch/mioverto/LOH_methods/Overton_etal/variants/final/Dual/all_BYxRM_hetTrain.RData"
# model_type <- "v"
# FE_rate <- 0.05
# SNPs_output <- "/tscc/lustre/ddn/scratch/mioverto/LOH_methods/Loeillet_etal/variants/final/Dual/all_BYxSK1_SNPs.RData"

SNPs_input <- as.character(input_args[1])
model_input <- as.character(input_args[2])
model_type <- as.character(input_args[3])
FE_rate <- as.numeric(input_args[4])
SNPs_output <- as.character(input_args[5])

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

cat("\nSNPs table imported with ", nrow(SNPs_merge_finalGT), " rows\n")


if(!any(c(grepl("Sum_DP", colnames(SNPs_merge_finalGT), 
        grepl("f_Alt", colnames(SNPs_merge_finalGT)))))) {
    SNPs_merge_finalGT <- SNPs_merge_finalGT %>% AddDPinfo()
}

load(file = model_input)

model_table <- model_list$model
varAlt_sumDP <- model_list$DP_table

cat("\nModel parameters imported. Beta = ", model_table$beta, "\n")

# End Section #################################################################
###############################################################################
#                                                                             #
#                                                                             #
###############################################################################
# Load in SNP data frame and model list #######################################

# Choose which model type to use
if(grepl("v", model_type, ignore.case = T)) {
    # Error beta parameter is used to predict variance as a function of Sum_DP:
    # var_Alt = 0.25 * Sum_DP + beta * Sum_DP^2
    # The alpha parameter of the beta-binomial distribution
    # (assuming BB_alpha = BB_beta) has an analytical solution
    # Given the BB_alpha calculated from the predicted variance, the upper and
    # lower quatiles for a beta-binomial at CDF < FE_rate/2 and 
    # (1 - CDF) < FE_rate/2. Function BB_rejectHets returns the 
    # set of calls that are het and outside the BB quantiles
    cat("\nFiltering based on predicted variance\n")
    SNPs_merge_finalGT$BB_reject <- SNPs_merge_finalGT %>% 
                                  select(GT, Alt_DP_final, Sum_DP_final) %>% 
                                  BB_rejectHets(., 
                                    error_beta = model_table$beta, 
                                    FE_rate = FE_rate)

    cat("\nFiltering complete\n")
} else if(grepl("d", model_type, ignore.case = T) |
            grepl("n", model_type, ignore.case = T) ) {
    cat("\nFiltering based on fixed Beta-binomial alpha\n")
    .alpha <- apply(varAlt_sumDP[, "BB_alpha"], MARGIN = 2, 
                    FUN = ifelse(grepl("d", model_type, ignore.case = T),
                                    median, mean))

    SNPs_merge_finalGT$BB_reject <- SNPs_merge_finalGT %>%
                                  select(GT, Alt_DP_final, Sum_DP_final) %>% 
                                  BB_rejectHets(., BB_alpha = .alpha,
                                  FE_rate = FE_rate)
} else {
    cat("Model type input not recognized")
}


output_RData <- grepl(".RData", SNPs_output)
if(output_RData) {
    save(SNPs_merge_finalGT, file = SNPs_output)
} else {
   write.table(SNPs_merge_finalGT, file = SNPs_output, sep = ",",
             col.names = T, row.names = F)
}

}
