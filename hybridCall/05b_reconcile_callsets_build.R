#!/usr/bin/env Rscript

###############################################################################
# This script loads in a VCF file and chain file for each parental callset.
# It parses the VCFs into data.frames, adjusts position indicies to a common 
# Type reference, and extracts SNPs. It then reconciles the two callsets at 
# each SNP position based on genotype likelihoods (GQ scores).
# The script requires functions from the 05a_reconcile_functions.R script.
# The expected VCF format is a multi-sample VCF output from GATK, and may 
# contain a combination of founders and/or end-point clones.
# The output is a single table of reconciled calls with column headers.
###############################################################################
###############################################################################
# Set parameters

# Minimum GQ for a reconciled call to be inlcuded
min_GQ_reconciled <- 50
# Minimum GQ to include a "winning" call when the other call is missing
min_GQ_pairNA <- 100
# Minimum difference in GQ values between "winning" and "losing" calls to retain
# the "winning" call
min_GQ_differential <- 30  

# Minimum read support for each allele
min_DP_per_allele <- 0

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
    make_option(c("-o", "--parentOne"), type = "character", default = NULL,
                action = "store", help = "Parent 1 VCF file name"),
    make_option(c("-t", "--parentTwo"), type = "character", default = NULL,
                action = "store", help = "Parent 2 VCF file name"),
    make_option(c("-c", "--parentOneChain"), type = "character",
                action = "store", default = NULL, 
                help = "Parent 1 chain file name"),
    make_option(c("-h", "--parentTwoChain"), type = "character",
                action = "store", default = NULL, 
                help = "Parent 2 chain file name"),
    make_option(c("-r", "--repeatBed"), type = "character", default = NULL,
                action = "store", help="Repeat .bed file name"),
    make_option(c("-F", "--TableFile"), type = "character", default = NULL, 
                action = "store", help = "output file (must be .RData or .csv)")
) 

# Parse arguements and create list of values
opt_parser <- OptionParser(option_list = option_list, add_help_option = FALSE)
input_args <- parse_args(opt_parser)

# Assign aguement list entries to separate variables
# P1_file <- as.character(input_args[1])
# P2_file <- as.character(input_args[2])
# P1_chainFile <- as.character(input_args[3])
# P2_chainFile <- as.character(input_args[4])
# repeats_bed <- as.character(input_args[5])
# output_table <- as.character(input_args[6])

samples <- "all"
P1_file <- paste0("/oasis/tscc/scratch/mioverto/LOH_methods/Overton_etal/variants/final/BY/", samples, "_BY.sif.vcf")
P2_file <- paste0("/oasis/tscc/scratch/mioverto/LOH_methods/Overton_etal/variants/final/RM/", samples, "_RM.sif.vcf")
P1_chainFile <- "/oasis/tscc/scratch/mioverto/LOH_methods/Overton_etal/references/final/BY/BYxS288C.chain"
P2_chainFile <- "/oasis/tscc/scratch/mioverto/LOH_methods/Overton_etal/references/final/RM/RMxS288C.chain"
repeats_bed <- "/oasis/tscc/scratch/mioverto/LOH_methods/Overton_etal/references/final/S288C/repeats/S288C_repeats.bed"
output_table <- paste0("/oasis/tscc/scratch/mioverto/LOH_methods/Overton_etal/variants/final/", samples, "_BYxRM_reconciled_SNPs.RData")

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
# Import chain files and chromosome boundaries for each parent

P1toRef_chainList <- ChainToDFlist(P1_chainFile)
P2toRef_chainList <- ChainToDFlist(P2_chainFile)

P1toRef_chainDf <- P1toRef_chainList$chain
P2toRef_chainDf <- P2toRef_chainList$chain

chrom_bounds_Ref  <- P1toRef_chainList$Ref_chrom
chrom_bounds_P1  <- P1toRef_chainList$Qry_chrom
chrom_bounds_P2  <- P2toRef_chainList$Qry_chrom

cat("Successfully created liftover and chromosome tables\n")

# End Section #################################################################
###############################################################################
#                                                                             #
#                                                                             #
###############################################################################
###############################################################################
# Import VCF files for each parent

# Generate master SNP dataframe from all clone data aligned to each reference
merge_cols <- c("CHROM", "POS", "POSi", "QUAL", "REF", "ALT", "GT", "Ref_DP", "Alt_DP",
                "GQ", "ID")

by_cols <- c("CHROM", "POS", "POSi", "ID")

# Load in Parent 1 vcf file
all_var_P1 <- HC_multiVcfParse_allVar(P1_file)

# Translate positions for P1 callset to Type positions
if(any(grepl("S288C", P1_chainFile))) {
all_var_P1 <- all_var_P1 %>%
  ConstructLiftover(.,  chain_df = P1toRef_chainDf, 
                    chrom_from = chrom_bounds_P1,
                    chrom_to = chrom_bounds_Ref)
}

# Load Parent 2 callset
all_var_P2 <- HC_multiVcfParse_allVar(P2_file)

# Translate positions for P2 callset to Type positions
all_var_P2 <- all_var_P2 %>%
    ConstructLiftover(.,  chain_df = P2toRef_chainDf, 
                      chrom_from = chrom_bounds_P2,
                      chrom_to = chrom_bounds_Ref)

# Get indicies of SNPs for each callset
iSNP_P1 <- nchar(all_var_P1$REF) == 1 & nchar(all_var_P1$ALT) == 1 
iSNP_P2 <- nchar(all_var_P2$REF) == 1 & nchar(all_var_P2$ALT) == 1

# End Section #################################################################
###############################################################################
#                                                                             #
#                                                                             #
###############################################################################
###############################################################################
# Merge two parental callsets based on sample ID and genomic positions (POSi)

SNPs_merge <- merge(all_var_P1[iSNP_P1, ], 
                    all_var_P2[iSNP_P2, ], 
                    by = by_cols, all = T, 
                    suffixes = c("_P1call", "_P2call"), 
                    sort = F) %>% arrange(ID, POSi)

sum_n_calls <- nrow(all_var_P1[iSNP_P1, ]) + nrow(all_var_P2[iSNP_P2, ])
mean_n_calls <- sum_n_calls / 2
merge_n_calls <- nrow(SNPs_merge)
merge_QC <- which.min(c(abs(merge_n_calls - sum_n_calls), 
                        abs(merge_n_calls - mean_n_calls)))

if(merge_QC == 1) {
    print("Merged callset length close to sum of individual callset legnths and 
           may be malformed")
}

rm(all_var_P1)
rm(all_var_P2)
rm(iSNP_P1)
rm(iSNP_P2)
gc()

# End Section #################################################################
###############################################################################
#                                                                             #
#                                                                             #
###############################################################################
# Correct allele order differences between callsets ###########################

SNPs_merge <- SwapCrossedAlleles(SNP_df = SNPs_merge,
                                 to_match = "P1call", to_swap = "P2call") %>%
                filter(!bad_alleles) %>%
                select(-bad_alleles)

# End Section #################################################################
###############################################################################
#                                                                             #
#                                                                             #
###############################################################################
# Remove sites that are in repeat regions #####################################

repeats_df <- ImportBed(repeats_bed)
SNPs_merge$repeats <- MarkRepeats(SNPs_merge, repeats_in = repeats_df)

cat("Fraction of calls in repeat regions", round(sum(SNPs_merge$repeats)/nrow(SNPs_merge), 2), "\n")

SNPs_merge <- SNPs_merge %>% filter(!repeats) %>% select(-repeats)

# End Section #################################################################
###############################################################################
#                                                                             #
#                                                                             #
###############################################################################
# Reconcile genotypes based on GQ values ######################################

SNPs_merge_finalGT <- SNPs_merge %>% 
                        GenotypeFromGQ(., P1_ID = "P1call", P2_ID = "P2call",
                            baseThrsh = min_GQ_reconciled, naThrsh = min_GQ_pairNA, 
                            diffThrsh = min_GQ_differential) %>%
                        filter(GT != "./.")

rm(SNPs_merge)

output_RData <- grepl(".RData", output_table)
if(output_RData) {
    save(all_var_P1, file = output_table)
} else {
   write.table(all_var_P1, file = output_table,
               sep = ",", col.names = T, row.names = F)
}
}