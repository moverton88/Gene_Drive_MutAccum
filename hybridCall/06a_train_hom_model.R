###############################################################################
# Even after the variant filtration steps included in GATK best practices,
# There are still homozygous calls that appear to be dubious given the
# excess in fraction of reads supporting an alternative allele. 
# However, we found that we could model the accummulation of erroneous Alt
# reads as a Poisson process. We do this by fitting a Poisson to k bins of Alt
# read counts such that x = 0:k encompasses at least 90% of the data. We then
# calculate the Alt read thresholds, above which we reject the call as too
# error-prone to yield a reliable call. 
###############################################################################
# Set parameters ##############################################################

# Minimum sample size for calcluating Alt read variance
min_n <- 100

# Max total depth to include in the variance vs depth model
# Total depth <= 200 includes ~98% of calls
model_max_DP <- 200

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

lib_dir <- "~/R/x86_64-pc-linux-gnu-library/3.6"
## Load packages
library(optparse, lib.loc = lib_dir)
library(stringr, lib.loc = lib_dir)
library(rlang, lib.loc = lib_dir)
library(emdbook, lib.loc = lib_dir)
library(vctrs, lib.loc = lib_dir)
library(cli, lib.loc = lib_dir)
library(tidyr, lib.loc = lib_dir)
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
    make_option(c("-h", "--homSitesIn"), type = "character", default = NULL,
                action = "store", help = "Homozygous sites Allele Depth table input"),
    make_option(c("-c", "--chainFile"), type = "character", default = NULL,
                action = "store", 
                help = "Chain file for lifting over positions to type reference"),
    make_option(c("-n", "--outputName"), type = "numeric", default = 0.9,
                action = "store", 
                help = "Name for output objects"),
    make_option(c("-o", "--outputDir"), type = "numeric", default = 0.5,
                action = "store", 
                help = "Output directory")
    make_option(c("-O", "--outputFormat"), type = "character", default = "R", 
                action = "store", help = "output format: R for .RData list object, C for two .csv tables")
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
# homAD_dir <- "/oasis/tscc/scratch/mioverto/LOH_methods/Overton_etal/variants/hom_AD"
# chain_file <- "/oasis/tscc/scratch/mioverto/LOH_methods/Overton_etal/references/final/BY/BYxS288C.chain"
# output_dir <- "/oasis/tscc/scratch/mioverto/LOH_methods/Overton_etal/variants/reconciled/"
# output_format <- "R"
homAD_dir <- as.character(input_args[1])
chain_file <- as.character(input_args[2])
output_name <- as.character(input_args[3])
output_dir <- as.character(input_args[4])
output_format <- as.character(input_args[5])

output_dir <- ifelse(grepl( "/$", output_dir), output_dir, paste0(output_dir, "/"))

# End Section #################################################################
###############################################################################
#                                                                             #
#                                                                             #
###############################################################################
# Import chain files and chromosome boundaries for each parent ################

P1toRef_chainList <- ChainToDFlist(chain_file)

P1toRef_chainDf <- P1toRef_chainList$chain

chrom_bounds_Ref  <- P1toRef_chainList$Ref_chrom
chrom_bounds_P1  <- P1toRef_chainList$Qry_chrom

cat("Successfully created liftover and chromosome tables\n")

# End Section #################################################################
###############################################################################
#                                                                             #
#                                                                             #
###############################################################################
# Import hom site allele depths ###############################################

hom_sites_AD <- data.frame()
homAD_dir_list <- list.files(homAD_dir)
for(i in homAD_dir_list) {
  # i <- anc_ID[1]
  homAD_file <- paste0(homAD_dir, "/", i)
  id_ADtable <- read.table(file = homAD_file, header = T, sep ="\t", 
                          col.names = c("CHROM", "POS", "GT", "AD", "GQ")) %>%
                separate(GT, into = c("REF", "ALT"), sep = "/",  extra = "drop") %>%
                separate(AD, into = c("Ref_DP", "Alt_DP"), sep = ",", convert = T, extra = "drop") %>%
                mutate(Sum_DP = Ref_DP + Alt_DP,
                        f_Alt = Alt_DP/Sum_DP,
                        ID = sample_ID)
  
  hom_sites_AD <- rbind(hom_sites_AD, id_ADtable)
}

hom_sites_AD <- hom_sites_AD %>% 
                  filter(Sum_DP > 1, GQ > 0) %>%
                  group_by(ID) %>%
                  distinct(POS, .keep_all = T) %>%
                  ungroup() %>%
                  mutate(CHROM = factor(CHROM),
                         ID = factor(ID))

hom_sites_AD <- hom_sites_AD %>% 
  mutate(CHROM = factor(chrom_bounds_Ref$CHROM[as.numeric(CHROM)])) %>%
  mutate(POSi = ConvertPosIndicies(., pos_col = "POS", chrom_col = "CHROM", 
                     index_out = "POSi", chrom_bound = chrom_bounds_P1,
                     add_chroms = F)) %>%
  ConstructLiftover(., chain_df = P1toRef_chain_df, 
                    pos_col = "Ref_POS", csDiff_col = "cs_Ref_diff")

dub_hom_posi <- hom_sites_AD %>% 
  filter(REF != ALT) %>% 
  distinct(POSi) %>% 
  pull(POSi)

hom_sites_AD <- hom_sites_AD %>% 
  filter(!POSi %in% dub_hom_posi)

rRNA_POSi_seq <- hom_sites_AD %>% distinct(POSi) %>% 
  filter(POSi >= rRNA_POSi[1], POSi <= rRNA_POSi[2]) %>% pull(POSi)

hom_sites_AD$repeats <- MarkRepeats(hom_sites_AD, repeats_in = repeats_bed)

hom_sites_AD <- hom_sites_AD %>% 
                  filter(!POSi %in% rRNA_POSi_seq &
                          !repeats) %>% 
                  select(-repeats)

# End Section #################################################################
###############################################################################
#                                                                             #
#                                                                             #
###############################################################################
# Filter to high-confidence homozygous sites ##################################

dub_hom_pos <- hom_sites_AD %>% 
                filter(REF != ALT) %>% 
                distinct(POS) %>% 
                pull(POS)

hom_sites_AD <- hom_sites_AD %>% 
                  filter(!POS %in% dub_hom_pos)

hom_sites_AD <- hom_sites_AD %>% 
  mutate(CHROM = factor(chrom_IDs$CHROM[as.numeric(CHROM)])) %>%
  mutate(POSi = ConvertPosIndicies(., pos_col = "POS", chrom_col = "CHROM", 
                     index_out = "POSi", chrom_bound = chrom_bounds_P1,
                     add_chroms = F)) %>%
  ConstructLiftover(., chain_df = P1toRefChainDf, pos_col = "Ref_POS", csDiff_col = "cs_Ref_diff")

hom_sites_AD$repeats <- MarkRepeats(hom_sites_AD, repeats_in = repeats_bed)

hom_sites_AD <- hom_sites_AD %>% 
                  filter(!repeats) %>% 
                  filter(!(POSi > rRNA_POSi[1] & POSi < rRNA_POSi[2])) %>% 
                  select(-repeats)

# End Section #################################################################
###############################################################################
#                                                                             #
#                                                                             #
###############################################################################
# Estimate Lambda from incremental distribution fraction ######################

mean_Lambdas_k2 <- hom_sites_AD %>% 
                    group_by(Sum_DP) %>% 
                    mutate(n = n()) %>% 
                    filter(n >= 100) %>%
                    summarize(L = LambdaFromsetKs(Alt_DP, k = 0:2))

# We limit Alt depth bin, k, to <= 2
# This is expected to include > 90% of the Alt depth distribution across 
# > 90% of total depth bins
Lambda_vs_Sum_DP_lm <- mean_Lambdas_k2 %>% 
                            lm(log(L) ~ log(Sum_DP), data = .)

model_coeffs <- summary(Lambda_vs_Sum_DP_lm)$coefficients
model_table <- data.frame(param = c("ln(alpha)", "beta"), model_coeffs, row.names = NULL)

colnames(model_table) <- c("param", "value", "StdErr", "t-value", "p-value")
model_table <- cbind(model_table, 
                     R_sq = summary(Lambda_vs_Sum_DP_lm)$adj.r.squared,
                     min_DP = min(mean_Lambdas_k2$Sum_DP),
                     max_DP = max(mean_Lambdas_k2$Sum_DP))


cat("Model output:\nAlpha = ", exp(model_table$value[1]) "Beta = ", model_table$value[2], 
    "\nR-squared = ", model_table$R_sq[1], "\n")

model_list <- list(hom_AltDP_df = hom_sites_AD, summary_table = model_table)

# End Section #################################################################
###############################################################################
#                                                                             #
#                                                                             #
###############################################################################
# Write hom sites table and model parameters ##################################

if(grepl("R", output_format, ignore.case = T)) {
    output_name <- paste0(output_name, "_homTrain.RData" )
    output_path <- paste0(output_dir, output_name)
    save(model_list, file = output_path)
} else {
    output_var_path <- paste0(output_dir, output_name, "_homAltDP.csv")
    output_model_path <- paste0(output_dir, output_name, "_homModel.csv")
    write.table(varAlt_sumDP, file = output_var_path, sep = ",",
            col.names = T, row.names = F)
    write.table(model_table, file = output_model_path, sep = ",",
            col.names = T, row.names = F)
}

}
