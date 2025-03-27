###############################################################################
# Even after the variant filtration steps included in GATK best practices,
# There are still heterozygous calls that appear to be dubious given the
# deviation in fraction of reads supporting the alternative allele from the 
# expected 0.5. In fact, if read generation follows a binomial process
# then we expect the variance in Alt read counts to be 0.25 * total reads.
# However, even if read generation is perfectly binomial, processing steps, 
# such as mapping and read filtering, could increase this variance. We would 
# like a model of these processes so that we can distinguish true from false
# heterozygous calls. 
###############################################################################
# Set parameters ##############################################################

# Min total depth to include in the variance vs depth model
# Depths under 6 produce BB alphas that result in thresholds
# < 1 or > Sum_DP
model_min_DP <- 6

# Max total depth to include in the variance vs depth model
# Total depth <= 200 includes ~98% of calls
# model_max_DP <- 200

# Minimum sample size for calcluating Alt read variance
# min_n <- 100

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
    make_option(c("-s", "--SNPsIn"), type = "character", default = NULL,
                action = "store", help = "SNPs table input"),
    make_option(c("-I", "--trainIDs"), type = "character", default = NULL,
                action = "store", 
                help = "Sample IDs to include in model training (string to grep or file with ID on each line)"),
    make_option(c("-f", "--fractionTrainIDs"), type = "numeric", default = 0.9,
                action = "store", 
                help = "Minimum fraction of training IDs with het call to include in training set"),
    make_option(c("-F", "--fractionNonTrainIDs"), type = "numeric", default = 0.5,
                action = "store", 
                help = "Minimum fraction of non-training IDs with het call to include in training set"),
    make_option(c("-p", "--inputParameters"), type = "character", default = 0.5,
                action = "store", 
                help = "Input parameters for the model (comma-sparated): max depth,"),
    make_option(c("-o", "--outputDir"), type = "character", default = NULL, 
                action = "store", help = "output directory"),
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

# Assign agruement list entries to separate variables

# SNPs_input <- "/oasis/tscc/scratch/mioverto/LOH_methods/Dutta_etal/variants/reconciled/H09_BAMxCPG_reconciled_SNPs.RData"
# train_ID <- "H09"
# fractionTrainIDs <- 0.9
# fractionNonTrainIDs <- 0.5
# input_params <- "4,1000,100"
# output_dir <- "/oasis/tscc/scratch/mioverto/LOH_methods/Dutta_etal/variants/reconciled/"
# output_format <- "R"

SNPs_input <- as.character(input_args[1])
train_ID <- as.character(input_args[2]) # Which clone IDs to include in model training (usually founder IDs) 
fractionTrainIDs <- as.numeric(input_args[3]) # Fraction of training clones that must have het call to call "known" het site
fractionNonTrainIDs <- as.numeric(input_args[4]) # Fraction of non-training clones that must have het call to call "known" het site
input_params <- as.character(input_args[5]) # Three parameters: min depth, max depth, and min calls in a depth bin for variance calculation
output_dir <- as.character(input_args[6])
output_format <- as.character(input_args[7])

input_params <- as.numeric(strsplit(input_params, ",", fixed = T)[[1]])
model_min_DP <- input_params[1]
model_max_DP <- input_params[2]
min_n <- input_params[3]

output_dir <- ifelse(grepl( "/$",output_dir), output_dir, paste0(output_dir, "/"))

# End Section #################################################################
###############################################################################
#                                                                             #
#                                                                             #
###############################################################################
# Load in SNP data frame and training IDs #####################################

input_RData <- grepl(".RData", SNPs_input)
if(input_RData) {
    load(file = SNPs_input)
} else {
   SNPs_merge_finalGT <- read.table(file = SNPs_input, header = T, sep = ",")
}

if(any(grepl("/", train_ID))) {
    train_IDs <- read.table(file = train_ID, stringsAsFactors = F) %>% pull(V1)
} else {
    train_IDs <- grep(train_ID, SNPs_merge_finalGT %>% distinct(ID) %>% pull(ID), 
                        value = T)
}

n_train_IDs <- fractionTrainIDs * length(train_IDs)
non_train_IDs <- SNPs_merge_finalGT %>% distinct(ID) %>% filter(!ID %in% train_IDs) %>% pull(ID)
n_non_train_IDs <- fractionNonTrainIDs * length(non_train_IDs)

# End Section #################################################################
###############################################################################
#                                                                             #
#                                                                             #
###############################################################################
# Get known heterozygous sites ################################################

if(!any(c(grepl("Sum_DP", colnames(SNPs_merge_finalGT), 
        grepl("f_Alt", colnames(SNPs_merge_finalGT)))))) {
    SNPs_merge_finalGT <- SNPs_merge_finalGT %>% AddDPinfo()
}

# Calculate site-wise genotype counts
sitewise_GTs <- SNPs_merge_finalGT %>%
  ungroup() %>%
  CalcGenoStats(., group = "all", anc_id = train_IDs)

# get known het sites given ancestral and end-point thresholds
known_het_POSi <- sitewise_GTs %>% 
                    filter(nHet_anc >= n_train_IDs, 
                           nHet_evo >= n_non_train_IDs) %>% 
                    pull(POSi)

cat("Number of known het sites = ", length(known_het_POSi), "\n")

# End Section #################################################################
###############################################################################
#                                                                             #
#                                                                             #
###############################################################################
# Calculate Alt read variance for each depth bin ##############################

# Calculate Alt read variance and sample size for each bin
# Must require >= 6 reads, since below this, variance is high and df is low
varAlt_sumDP <- SNPs_merge_finalGT %>% 
  filter(ID %in% train_IDs, GT == "0/1", 
            Sum_DP_final >= model_min_DP, POSi %in% known_het_POSi) %>%
  select(POSi, GT, GQ, Alt_DP_final, Sum_DP_final, f_Alt) %>%
  group_by(Sum_DP_final) %>%
  summarize(n = n(),
            Sum_DP_mean = mean(Sum_DP_final),
            mean_Alt_DP = mean(Alt_DP_final),
            var_Alt_DP = var(Alt_DP_final),
            var_f_Alt = var(f_Alt))

# Get total depths with sample sizes over a threshold
valid_DP_values <- varAlt_sumDP %>% 
                    filter(n >= min_n, 
                           Sum_DP_mean >= model_min_DP, 
                           Sum_DP_mean <= model_max_DP) %>% 
                    pull(Sum_DP_mean)
cat("Number of depth bins for model = ", length(valid_DP_values), "\n")
# End Section #################################################################
###############################################################################
#                                                                             #
#                                                                             #
###############################################################################
# Fixed binomial model with quadratic error term ##############################

varAlt_sumDP <- varAlt_sumDP %>% 
                mutate(var_Alt_resid = var_Alt_DP - 0.25 * Sum_DP_mean)

if(length(valid_DP_values) > 10) {
    var_lm_fix_Binom <- varAlt_sumDP %>% 
                            filter(Sum_DP_mean %in% valid_DP_values) %>%
                            mutate(BB_nSqMinusn = Sum_DP_mean^2 - Sum_DP_mean) %>%
                            lm(var_Alt_resid ~ 0 + BB_nSqMinusn, data = .)

    var_model <- c(0.25, var_lm_fix_Binom$coefficients)

    varAlt_sumDP <- varAlt_sumDP %>%
        mutate(pred_Alt_var = var_model[1] * Sum_DP_mean + 
                var_model[2] * Sum_DP_mean ^ 2,
                pred_BB_alpha = BB_alphaFromVar(pred_Alt_var, Sum_DP_mean),
                obs_BB_alpha = BB_alphaFromVar(var_Alt_DP, Sum_DP_mean))

    model_table <- summary(var_lm_fix_Binom)$coefficients %>% 
                    data.frame(row.names = 1)
    
    colnames(model_table) <- c("beta", "StdErr", "t-value", "p-value")
    model_table <- cbind(model_table, 
                     R_sq = summary(var_lm_fix_Binom)$adj.r.squared,
                     min_DP = min(valid_DP_values),
                     max_DP = max(valid_DP_values))
} else {
    var_model <- c(0.25, NA)

    varAlt_sumDP <- varAlt_sumDP %>%
        mutate(pred_Alt_var = NA,
                pred_BB_alpha = NA,
                obs_BB_alpha = NA)

    model_table <- data.frame(beta = NA)
}

cat("Model output:\nBeta = ", model_table$beta, 
    "\nR-squared = ", model_table$R_sq, "\n")

model_list <- list(var_Alt_table = varAlt_sumDP, summary_table = model_table)

# End Section #################################################################
###############################################################################
#                                                                             #
#                                                                             #
###############################################################################
# Write SNPs table and model parameters #######################################

output_name <- sub(".", "_", strsplit(sub("_", ".", basename(SNPs_input)), "_")[[1]][1], fixed = T)

if(grepl("R", output_format, ignore.case = T)) {
    output_name <- paste0(output_name, "_hetTrain.RData" )
    output_path <- paste0(output_dir, output_name)
    save(model_list, file = output_path)
} else {
    output_var_path <- paste0(output_dir, output_name, "_hetAltVar.csv")
    output_model_path <- paste0(output_dir, output_name, "_hetModel.csv")
    write.table(varAlt_sumDP, file = output_var_path, sep = ",",
            col.names = T, row.names = F)
    write.table(model_table, file = output_model_path, sep = ",",
            col.names = T, row.names = F)
}


}
