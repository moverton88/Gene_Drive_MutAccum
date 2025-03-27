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
min_GQ_reconciled <- 30
# Minimum GQ to include a "winning" call when the other call is missing
min_GQ_pairNA <- 100
# Minimum difference in GQ values between "winning" and "losing" calls to retain
# the "winning" call
min_GQ_differential <- 30

# Minimum read support for each allele
min_DP_per_allele <- 0

# Minimum read support for each allele
min_QUAL <- 1000

# Chr XII rRNA array genomic coordinates
# 7245029 + c(447000, 490000)
rRNA_POSi <- c(7692029, 7735029)

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
    make_option(c("-S", "--SNPsOut"), type = "character", default = NULL, 
                action = "store", help = "output SNP table file (must be .RData or .csv)")
) 

# Parse arguements and create list of values
opt_parser <- OptionParser(option_list = option_list, add_help_option = FALSE)
input_args <- parse_args(opt_parser)

# Assign aguement list entries to separate variables
# P1_file <- "/tscc/lustre/ddn/scratch/mioverto/LOH_methods/Pankajam_etal/variants/final/S288C/L_S288C.filter.vcf"
# P2_file <- "/tscc/lustre/ddn/scratch/mioverto/LOH_methods/Pankajam_etal/variants/final/RM/L_RM.filter.vcf"
# P1_chainFile <- ""
# P2_chainFile <- "/tscc/lustre/ddn/scratch/mioverto/LOH_methods/Pankajam_etal/references/final/RM/RMxS288C.chain"
# repeats_bed <- "/tscc/lustre/ddn/scratch/mioverto/LOH_methods/Pankajam_etal/references/final/S288C/repeats/S288C_repeats.bed"
# output_table <- "/tscc/lustre/ddn/scratch/mioverto/LOH_methods/Pankajam_etal/variants/final/S288C/L_S288CxRM_reconciled_SNPs.RData"
# dual_file_old <- "/tscc/lustre/ddn/scratch/mioverto/LOH_methods/Pankajam_etal/variants/final/Dual/L_S288CxRM_SNPs.RData"
# dual_file <- "/tscc/lustre/ddn/scratch/mioverto/LOH_methods/Pankajam_etal/variants/final/Dual/L_S288CxRM_reconciled_SNPs.RData"

P1_file <- as.character(input_args[1])
P2_file <- as.character(input_args[2])
P1_chainFile <- as.character(input_args[3])
P2_chainFile <- as.character(input_args[4])
repeats_bed <- as.character(input_args[5])
output_table <- as.character(input_args[6])

# Determine whether arguements are valid and stop script if not
args_invalid <- any(sapply(input_args, function(x) is.null(x)))
if(args_invalid) {
    print("One or more options not provided")
} else {

P1_filename <- str_split(P1_file, "") %>% unlist()
i_undscr <- grep("_", P1_filename)
i_dot <- grep('\\.', P1_filename)
P1 <- P1_filename[(i_undscr[length(i_undscr)] + 1):(i_dot[1] - 1)] %>% paste(collapse = "")
P2_filename <- str_split(P2_file, "") %>% unlist()
i_undscr <- grep("_", P2_filename)
i_dot <- grep('\\.', P2_filename)
P2 <- P2_filename[(i_undscr[length(i_undscr)] + 1):(i_dot[1] - 1)] %>% paste(collapse = "")

# End Section #################################################################
###############################################################################
#                                                                             #
#                                                                             #
###############################################################################
###############################################################################
# Import chain files and chromosome boundaries for each parent
if(any(grepl("xS288C", P1_chainFile))) {
    print("P1 is not strain S288C, attempting to open P1 chain file")
    P1toRef_chainList <- ChainToDFlist(P1_chainFile)
    chrom_bounds_P1  <- P1toRef_chainList$Qry_chrom
    P1toRef_chainDf <- P1toRef_chainList$chain
}

P2toRef_chainList <- ChainToDFlist(P2_chainFile)
P2toRef_chainDf <- P2toRef_chainList$chain
chrom_bounds_P2  <- P2toRef_chainList$Qry_chrom
chrom_bounds_Ref  <- P2toRef_chainList$Ref_chrom

cat("Successfully created liftover and chromosome tables\n")

# End Section #################################################################
###############################################################################
#                                                                             #
#                                                                             #
###############################################################################
###############################################################################
# Import VCF files for each parent

# Generate master SNP dataframe from all clone data aligned to each reference
merge_cols <- c("CHROM", "POS", "POSi", "QUAL", "REF", 
                "ALT", "GT", "Ref_DP", "Alt_DP", "Sum_DP",
                "GQ", "ID")
by_cols <- c("CHROM", "POS", "POSi", "ID")

# Load in Parent 1 vcf file
all_var_P1 <- HC_multiVcfParse_allVar(P1_file)
all_var_P1$ID <- gsub(paste0("_", P1), "", all_var_P1$ID)

cat("Parent 1 VCF read in\n")
# Translate positions for P1 callset to Type positions
if(any(grepl("xS288C", P1_chainFile))) {
all_var_P1_lo <- all_var_P1 %>% 
  ConstructLiftover(.,  chain_df = P1toRef_chainDf, 
                    chrom_from = chrom_bounds_P1,
                    chrom_to = chrom_bounds_Ref)
    cat("Parent 1 positions lifted over to S288C\n")
}

# Load Parent 2 callset
all_var_P2 <- HC_multiVcfParse_allVar(P2_file)
all_var_P2$ID <- gsub(paste0("_", P2), "", all_var_P2$ID)
cat("Parent 2 VCF read in")
# Translate positions for P2 callset to Type positions
all_var_P2_lo <- all_var_P2 %>%
    ConstructLiftover(.,  chain_df = P2toRef_chainDf, 
                      chrom_from = chrom_bounds_P2,
                      chrom_to = chrom_bounds_Ref)
cat("Parent 2 positions lifted over to S288C")

P1_1 <- all_var_P1_lo %>% filter(ID == "wt_00")
P2_1 <- all_var_P1_lo %>% filter(ID == "wt_00")

# Get indicies of SNPs for each callset
iSNP_P1 <- nchar(all_var_P1_lo$REF) == 1 & nchar(all_var_P1_lo$ALT) == 1 
iSNP_P2 <- nchar(all_var_P2_lo$REF) == 1 & nchar(all_var_P2_lo$ALT) == 1

# End Section #################################################################
###############################################################################
#                                                                             #
#                                                                             #
###############################################################################
###############################################################################
# Merge two parental callsets based on sample ID and genomic positions (POSi)

SNPs_merge <- merge(all_var_P1_lo[iSNP_P1, ], 
                    all_var_P2_lo[iSNP_P2, ], 
                    by = by_cols, all = T, 
                    suffixes = c("_P1call", "_P2call"), 
                    sort = F) %>% arrange(ID, POSi) %>%
                    mutate(ID = factor(ID))

n_P1calls <- sum(!is.na(SNPs_merge$GT_P1call))
n_P2calls <- sum(!is.na(SNPs_merge$GT_P2call))
n_discord <- sum(is.na(SNPs_merge$GT_P1call) | is.na(SNPs_merge$GT_P2call))
sum_n_calls <- (n_P1calls + n_P2calls)
mean_n_calls <- sum_n_calls / 2
merge_n_calls <- nrow(SNPs_merge)
merge_QC <- which.min(c(abs(merge_n_calls - sum_n_calls), 
                        abs(merge_n_calls - mean_n_calls)))

if(merge_QC == 1) {
    cat("Merged callset length close to sum of individual callset legnths and 
           may be malformed")
}

reconcile_stats <- data.frame(n_P1call = n_P1calls,
                              n_P2call = n_P2calls,
                              n_merge = merge_n_calls)

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
                                 to_match = "P1call", to_swap = "P2call")

# End Section #################################################################
###############################################################################
#                                                                             #
#                                                                             #
###############################################################################
# Remove sites that are in repeat regions #####################################

repeats_df <- ImportBed(repeats_bed)
SNPs_merge$repeats <- MarkRepeats(SNPs_merge, repeats_in = repeats_df)

cat("Fraction of calls in repeat regions", round(sum(SNPs_merge$repeats)/nrow(SNPs_merge), 2), "\n")

SNPs_merge <- SNPs_merge %>% 
                filter(!repeats) %>%
                filter(!(POSi > rRNA_POSi[1] & POSi < rRNA_POSi[2])) %>% 
                select(-repeats)

# End Section #################################################################
###############################################################################
#                                                                             #
#                                                                             #
###############################################################################
# Reconcile genotypes based on GQ values ######################################

SNPs_merge_finalGT <- SNPs_merge %>% 
                        GenotypeFromGQ(., P1_ID = "P1call", P2_ID = "P2call",
                            baseThrsh = min_GQ_reconciled, naThrsh = min_GQ_pairNA, diffThrsh = min_GQ_differential) %>%
                        filter(GT != "./.") %>%
                        mutate(GT = factor(GT, levels = c("0/0", "0/1", "1/1"))) %>%
                        AddDPinfo()

rm(SNPs_merge)

# Remove all calls with insufficient support
# Het - at least three reads for each allele
# Hom - at least three reads supporting allele
# Fill in depth information
if(min_DP_per_allele > 0) {
    SNPs_merge_finalGT <- SNPs_merge_finalGT %>%
        filter(!(GT == "0/1" & (Ref_DP_final <= min_DP_per_allele |
                                Alt_DP_final <= min_DP_per_allele)),
            !(GT == "1/1" & Alt_DP_final <= min_DP_per_allele),
            !(GT == "0/0" & Ref_DP_final <= min_DP_per_allele))
}

if(min_QUAL > 0) {
    SNPs_merge_finalGT <- SNPs_merge_finalGT %>%
        filter(QUAL_P1call >= min_QUAL | QUAL_P2call >= min_QUAL)
}


reconcile_stats <- reconcile_stats %>%
                    mutate(n_reconciled = nrow(SNPs_merge_finalGT),
                            p_P1call_rec = n_P1calls/n_reconciled,
                            p_P2call_rec = n_P2calls/n_reconciled)


output_RData <- grepl(".RData", output_table)
if(output_RData) {
    save(SNPs_merge_finalGT, file = output_table)
} else {
   write.table(SNPs_merge_finalGT, file = output_table, sep = ",",
             col.names = T, row.names = F)
}

reconcile_path <- paste0(file.path(output_table), "reconcile_stats.csv")
write.table(reconcile_stats, file = , sep = ",",
             col.names = T, row.names = F)
}

# End Section #################################################################
###############################################################################
#                                                                             #
#                                                                             #
###############################################################################
# End Script ##################################################################
sub_var_P1 <- all_var_P1 %>% filter(ID == "wt_12")
sub_var_P2 <- all_var_P2 %>% filter(ID == "wt_12")

output_table <- "/tscc/lustre/ddn/scratch/mioverto/LOH_methods/Loeillet_etal/variants/final/Dual/sub_var_P1.RData"
output_table <- "/tscc/lustre/ddn/scratch/mioverto/LOH_methods/Loeillet_etal/variants/final/Dual/sub_var_P2.RData"
save(SNPs_merge_finalGT, file = output_table)
}