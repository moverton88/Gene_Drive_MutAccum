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

type_strain <- "S288C"
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
    make_option(c("-c", "--parentOneChain"), type = "character",
                action = "store", default = NULL, 
                help = "Parent 1 chain file name"),
    make_option(c("-r", "--repeatBed"), type = "character", default = NULL,
                action = "store", help="Repeat .bed file name"),
    make_option(c("-S", "--SNPsOut"), type = "character", default = NULL, 
                action = "store", help = "output SNP table file (must be .RData or .csv)")
) 

# Parse arguements and create list of values
opt_parser <- OptionParser(option_list = option_list, add_help_option = FALSE)
input_args <- parse_args(opt_parser)

# Assign aguement list entries to separate variables
# P1_file <- "/tscc/lustre/ddn/scratch/mioverto/LOH_methods/Overton_etal/variants/final/S288C/all_S288C.filter.vcf"
# P1_chainFile <- "/tscc/lustre/ddn/scratch/mioverto/LOH_methods/Overton_etal/references/final/S288C/S288C_chrom_lengths.tsv"
# repeats_bed <- "/tscc/lustre/ddn/scratch/mioverto/LOH_methods/Overton_etal/references/final/S288C/repeats/S288C_repeats.bed"
# output_table <- "/tscc/lustre/ddn/scratch/mioverto/LOH_methods/Overton_etal/variants/final/S288C/all_S288C_singlRef_SNPs.RData"

P1_file <- as.character(input_args[1])
P1_chainFile <- as.character(input_args[2])
repeats_bed <- as.character(input_args[3])
output_table <- as.character(input_args[4])

# Determine whether arguements are valid and stop script if not
args_invalid <- any(sapply(input_args, function(x) is.null(x)))
if(args_invalid) {
    print("One or more options not provided")
} else {

P1_filename <- str_split(P1_file, "") %>% unlist()
i_undscr <- grep("_", P1_filename)
i_dot <- grep('\\.', P1_filename)
P1 <- P1_filename[(i_undscr[length(i_undscr)] + 1):(i_dot[1] - 1)] %>% paste(collapse = "")

# End Section #################################################################
###############################################################################
#                                                                             #
#                                                                             #
###############################################################################
###############################################################################
# Import chain files and chromosome boundaries for each parent
if(any(grepl("chain", P1_chainFile, ignore.case = T))) {
    cat("P1 is not strain S288C, attempting to open P1 chain file")
    P1toRef_chainList <- ChainToDFlist(P1_chainFile)
    chrom_bounds_Ref  <- P1toRef_chainList$Ref_chrom
    chrom_bounds_P1  <- P1toRef_chainList$Qry_chrom
    P1toRef_chainDf <- P1toRef_chainList$chain
} else if(any(grepl("length", P1_chainFile, ignore.case = T))){
    chrom_bounds_Ref  <- read.table(file = P1_chainFile, sep = "\t", col.names = c("CHROM_name", "length")) %>%
                            head(16) %>%
                            mutate(CHROM = str_pad(row_number(), 2, pad = "0")) %>%
                            select(CHROM, length)
} else {
    cat("Chain file/chromosome table option not recognized")
}

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
if(!any(grepl(paste0("_", type_strain), P1_file))) {
all_var_P1 <- all_var_P1 %>%
  ConstructLiftover(.,  chain_df = P1toRef_chainDf, 
                    chrom_from = chrom_bounds_P1,
                    chrom_to = chrom_bounds_Ref)
    cat("Parent 1 positions lifted over to S288C\n")
}


# Get indicies of SNPs for each callset
iSNP_P1 <- nchar(all_var_P1$REF) == 1 & nchar(all_var_P1$ALT) == 1 

# End Section #################################################################
###############################################################################
#                                                                             #
#                                                                             #
###############################################################################
###############################################################################
# Merge two parental callsets based on sample ID and genomic positions (POSi)

SNPs_merge <- all_var_P1 %>%
                filter(iSNP_P1) %>% 
                rename_with(.fn = ~ paste0(., "_final"), 
                            .cols = any_of(c("Ref_DP", "Alt_DP", "Sum_DP"))) %>% 
                            arrange(ID, POSi) %>%
                rename(QUAL_P1call = QUAL) %>%
                mutate(ID = factor(ID))



rm(all_var_P1)
rm(iSNP_P1)
gc()

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
                        filter(GQ >= min_GQ_reconciled, GT %in% c("0/0", "0/1", "1/1")) %>%
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
        filter(QUAL_P1call >= min_QUAL)
}

output_RData <- grepl(".RData", output_table)
if(output_RData) {
    save(SNPs_merge_finalGT, file = output_table)
} else {
   write.table(SNPs_merge_finalGT, file = output_table, sep = ",",
             col.names = T, row.names = F)
}
}

# End Section #################################################################
###############################################################################
#                                                                             #
#                                                                             #
###############################################################################
# End Script ##################################################################