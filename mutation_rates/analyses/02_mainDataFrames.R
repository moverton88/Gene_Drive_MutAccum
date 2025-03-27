###############################################################################
# Uses objects from 01_initialize_analyses and functions from 
# import_functions.R and analysis_functions.R to load in VCF files and 
# convert to data.frames.
# The expected format is a multi-sample VCF that contains all founders (anc) 
# and all end-point clones (evo)

# Depth stats from GATK CollectWgsMetrics for each clone and Reference
###############################################################################
# Import and bind depth metrics for all clones and each reference
BY_depths <- import_depths(depthDir = paste0(dataDir, "depth_metrics/"), callSet = "BY")
RM_depths <- import_depths(depthDir = paste0(dataDir, "depth_metrics/"), callSet = "RM")

# Get colnames to operate on
ID_cols <- c("ID", "Tx")
value_cols <- c("n_sites", "n_low_cover", 
                "n_valid", "MEDIAN_COVERAGE", "PCT_EXC_TOTAL")

# Get mean values between references and bind into one dataframe
depths_int <- do.call(rbind, list(BY_depths[, c(ID_cols, value_cols)], 
                                  RM_depths[, c(ID_cols, value_cols)]))
seq_depths_mean <- aggregate(depths_int[-c(1, 2)], 
                             list(depths_int$ID, depths_int$Tx), mean)
colnames(seq_depths_mean)[1:2] <- ID_cols
seq_depths_mean[, value_cols[-5]] <- lapply(seq_depths_mean[, value_cols[-5]], 
                                            round)
seq_depths_mean <- seq_depths_mean %>% mutate(f_valid = n_valid/g_length)
seq_depths_mean$Tx_ID <- seq_depths_mean$Tx %>% Recode_Tx_ID()
save(seq_depths_mean, file = depth_filename)

rm(BY_depths)
rm(RM_depths)

# Only retain clones with >= 70% coverage
low_depth_ID <- seq_depths_mean %>% 
  filter(n_valid/g_length < 0.7) %>%
  distinct(ID) %>% pull(ID)

# low_depth_ID <- c("F_A01", "F_D00", "F_E10", "F_E13", "F_G02",
#                   "H_A03", "H_A09", "H_E08", "H_E10", "H_F00", "H_H00",
#                   "N_A01", "N_A06", "N_A07", "N_A08", "N_C03")

# Calculate mean coverage
mean_cover <- seq_depths_mean %>% 
  filter(f_valid > 0.7) %>% 
  pull(n_valid) %>% mean()

# mean_cover <- 10446075

# Generate master SNP dataframe from all clone data aligned to each reference
###############################################################################

# Load in all vcf files and filter out low depth sites
all_var_BY <- HC_multiVcfParse_allVar(BY_file, all_Alts = T) %>% 
  filter(Sum_DP >= 6)

all_var_RM <- HC_multiVcfParse_allVar(RM_file, all_Alts = T) %>% 
  filter(Sum_DP >= 6)

# Translate positions for RM callset to BY positions
all_var_RM <- all_var_RM %>% 
  RMxBY_liftover(.,  chain_df = BYtoRMchainDf, lift_from = "RM", lift_to = "BY")

save(all_var_BY, file = BY_raw_filename)
save(all_var_RM, file = RM_raw_filename)

# Index putative biallelic SNPs by finding sites where the first two alleles are 
# single characters and the third allele is na 

iSNP_BY <- nchar(all_var_BY$REF) == 1 & nchar(all_var_BY$ALT1) == 1 
iSNP_RM <- nchar(all_var_RM$REF) == 1 & nchar(all_var_RM$ALT1) == 1

iBi_BY <- is.na(all_var_BY$ALT2)
iBi_RM <- is.na(all_var_RM$ALT2)

# Merge BY and RM dataframes by position and clone
merge_cols <- c("CHROM", "POS", "POSi", "QUAL", "REF", "ALT1", "GT", "Ref_DP", "Alt1_DP",
                "GQ", "Tx", "Line", "Rep", "ID")

by_cols <- c("CHROM", "POS", "POSi", "Tx", "Line", "Rep", "ID")

SNPs_merge <- merge(all_var_BY[iSNP_BY & iBi_BY, merge_cols], 
                    all_var_RM[iSNP_RM & iBi_RM, merge_cols], 
                    by = by_cols, all = T, 
                    suffixes = c("_BYcall", "_RMcall"), 
                    sort = F) %>% arrange(ID, POSi)

colnames(SNPs_merge)[(grep("ALT1", colnames(SNPs_merge)))] <- c("ALT_BYcall", "ALT_RMcall")
colnames(SNPs_merge)[(grep("Alt1", colnames(SNPs_merge)))] <- c("Alt_DP_BYcall", "Alt_DP_RMcall")

SNPs_merge <- SwapCrossedAlleles(SNP_df = SNPs_merge, annotated_POSi = RMxBY_comp_SNPs_POSi, 
                                 to_match = "BYcall", to_swap = "RMcall")

colnames(SNPs_merge) <- gsub("GenQual", "GQ", colnames(SNPs_merge))
SNPs_merge$Tx_name <- SNPs_merge$Tx %>% Recode_Tx()
SNPs_merge$Tx_ID <- SNPs_merge$Tx %>% Recode_Tx_ID()

rm(all_var_BY)
rm(all_var_RM)

# Remove sites that are in the rRNA repeat region and known Ty elements
###############################################################################

rRNA_POSi_seq <- SNPs_merge %>% distinct(POSi) %>% 
  filter(POSi >= rRNA_POSi$POSi[1], POSi <= rRNA_POSi$POSi[2]) %>% pull(POSi)
SNPs_merge <- SNPs_merge %>% filter(!POSi %in% rRNA_POSi_seq)

Ty_POSi_seq <- SNPs_merge %>% distinct(POSi) %>% 
  filter(POSi >= Ty_POSi$POSi[1], POSi <= Ty_POSi$POSi[2]) %>% pull(POSi)
SNPs_merge <- SNPs_merge %>% filter(!POSi %in% Ty_POSi_seq)

SNPs_merge$repeats <- MarkRepeats(SNPs_merge, repeats_in = repeats_bed)
SNPs_merge <- SNPs_merge %>% filter(!repeats) %>% select(-repeats)

## Write table of merged VCFs ##
# save(SNPs_merge, file = SNPs_merge_raw_filename)
# load(file = SNPs_merge_raw_filename)

# Number of sites initially sampled
n_BY_sites <- SNPs_merge %>% 
  filter(!is.na(GT_BYcall), !ID %in% low_depth_ID) %>% 
  distinct(POSi) %>% nrow()
n_RM_sites <- SNPs_merge %>% 
  filter(!is.na(GT_RMcall), !ID %in% low_depth_ID) %>% 
  distinct(POSi) %>% nrow()
n_sites_mean <- round(mean(c(n_BY_sites, n_RM_sites)))

n_sites_ID_BY <- SNPs_merge %>% 
  filter(!is.na(GT_BYcall), !ID %in% low_depth_ID) %>% 
  count(Tx, ID)
n_sites_ID_RM <- SNPs_merge %>% 
  filter(!is.na(GT_RMcall), !ID %in% low_depth_ID) %>% 
  count(Tx, ID)
n_sites_ID <- merge(n_sites_ID_BY, n_sites_ID_RM, 
                    by = c("Tx", "ID"), suffixes = c("_BY", "_RM")) %>%
  mutate(n_mean = round((n_BY + n_RM) / 2))
n_sites_ID$Tx_ID <- n_sites_ID$Tx %>% Recode_Tx_ID()

# some samples were poorly sequenced or lacked coverage across the entire genome 
# as well, DNA extraction failed for two WT founders 
bad_seq <- c(low_depth_ID, "N_F00", "N_H00")

# Preliminary LOH analysis identified clones with contamination
# Contiguous LOH tracts were interrupted by sporadic Het calls 
# with high allele imbalance favoring the allele of the LOH
contaminated <- c("F_D04", "N_B08", "N_D08", "N_D10", "N_E09", 
                  "N_E10", "N_F04", "N_G03", "H_D05", "H_B12", 
                  "H_E03", "H_G05", "H_G07")

# Some clones have whole genome aneuploidies and were removed from 
# LOH analysis
ID_allele_bias <- SNPs_merge %>% 
  filter(!ID %in% low_depth_ID, GT_BYcall == "0/1", GT_RMcall == "0/1") %>% 
  group_by(Line, ID) %>% 
  summarize(bias_Ref = sum(Ref_DP_BYcall + Ref_DP_RMcall)/
              sum(Ref_DP_BYcall + Alt_DP_BYcall + 
                    Ref_DP_RMcall + Alt_DP_RMcall))

ID_allele_bias %>% 
  ggplot() + geom_histogram(aes(x = bias_Ref), binwidth = 0.01) + 
  xlim(0, 1)

ID_allele_bias %>% 
  filter(bias_Ref < 0.45 | bias_Ref > 0.55) %>% 
  View()

whGnm_aneu_line <- c("F_B", "F_E")
whGnm_aneu_ID <- c("N_B12")

# Does the heterozygous call allele bias due to triploidy interfere with
# our ability to filter out dubious heterozygous calls?
# High-confidence het markers have about the same allele bias dist
# as unknown ones. Since we generate the confidence intervals with
# empirical distributions, this should work with bias due to 
# triploidy. 
SNPs_merge %>% 
  filter(Line %in% c("F_B", "F_E"), GT_RMcall == "0/1") %>% 
  ggplot() + 
  geom_histogram(aes(x = Ref_DP_RMcall/(Ref_DP_RMcall + Alt_DP_RMcall))) +
  scale_y_log10()

# These founder groups are processed in parallel with the script
# triploid_mainDataFrames.R


# Get the set of founders that did not sequence well for 
# later genotype imputation
anc_lines <- SNPs_merge %>% 
  filter(!(ID %in% c(contaminated, bad_seq, whGnm_aneu_ID))) %>% 
  filter(Rep == "00") %>% distinct(Line) %>% pull(Line)

all_lines <- SNPs_merge %>% distinct(Line)
noAncestor <- all_lines %>% 
  filter(!Line %in% anc_lines) %>% 
  pull(Line) %>% as.character()

# noAncestor <- c("F_D", "H_F", "H_H", "N_F", "N_H")

# Remove clones with confounding factors and perform final genotyping
SNPs_merge_finalGT <- SNPs_merge %>% 
  filter(!(ID %in% c(contaminated, bad_seq, whGnm_aneu_ID) | 
             Line %in% whGnm_aneu_line)) %>%
  GenotypeFromGQ(., baseThrsh = 50, naThrsh = 100, diffThrsh = 30) %>%
  filter(finalGT != "./.", finalGQ >= 50)

# with GQ > 50, nrow(SNPs_merge_finalGT) = 14843616 * 0.6372622 = 9459275
# with GQ > 30, nrow(SNPs_merge_finalGT) = 9844340 (0.6630044)
# with na GQ > 90, nrow(SNPs_merge_finalGT) = 10506275 (0.7077975)

# SNPs_merge object is no longer needed
# rm(SNPs_merge)
# load(file = SNPs_merge_raw_filename)

# Fill in information, standardize column names, drop levels, 
# find annotated SNPs from RMxBY tables
SNPs_merge_finalGT <- SNPs_merge_finalGT %>% FormatSNPsMerge()

# Remove all calls with insufficient support
# Het - at least three reads for each allele
# Hom - at least three reads supporting allele
SNPs_merge_finalGT <- SNPs_merge_finalGT %>%
  filter(!(GT == "0/1" & (Ref_DP_final <= 3 | Alt_DP_final <= 3)),
         !(GT == "1/1" & Alt_DP_final <= 3),
         !(GT == "0/0" & Ref_DP_final <= 3)) 

# Mark heterozygous calls that have excessive allele imbalance
SNPs_merge_finalGT$Cut <- MarkDubiousHets(SNPs_merge_finalGT, .alpha = 0.05)

# Impute genotypes of missing ancestors from end-point clone data
SNPs_merge_finalGT <- SNPs_merge_finalGT %>% 
  InferMissingFounders(missing_founders = noAncestor)
SNPs_merge_finalGT <- SNPs_merge_finalGT %>% FormatSNPsMerge()

# Are there any duplicated positions in any clones?
SNPs_merge_finalGT %>%
  count(ID, POSi) %>% 
  filter(n > 1) %>% nrow()

# load(file = SNPs_merge_filename)
# Count each genotype for each site across founder clones and 
# end-point clones independently
sitewise_GTs <- SNPs_merge_finalGT %>% 
  filter(!Cut) %>%
  site_genotype_stats(group = "all") %>%
  mutate(nHom_anc = nRef_anc + nRef_anc)

# Get IDs of all end-point clones with SNP data
all_IDs <- SNPs_merge_finalGT %>% 
  distinct(ID, .keep_all = T) %>% 
  pull(ID) %>% as.character()

evo_IDs <- SNPs_merge_finalGT %>% 
  distinct(ID, .keep_all = T) %>% 
  filter(Rep != "00") %>% pull(ID) %>% as.character()

n_founders <- SNPs_merge_finalGT %>% 
  filter(Rep == "00") %>%
  distinct(ID, .keep_all = T) %>% count(Tx_name)

n_clones <- SNPs_merge_finalGT %>% 
  filter(Rep != "00") %>% 
  pull(ID) %>% unique() %>% length()

n_clones_LOH_xTx <- SNPs_merge_finalGT %>% 
  distinct(ID, .keep_all = T) %>%
  filter(Rep != "00") %>%
  count(Tx_ID)

n_clones_LOH_xTx$Tx_ID <- Recode_Tx_ID(n_clones_LOH_xTx$Tx_ID, tx_type = "Tx_ID")
n_clones_LOH_xTx <- n_clones_LOH_xTx %>% arrange(Tx_ID)

n_all_xTx <- SNPs_merge_finalGT %>% 
  distinct(ID, .keep_all = T) %>% count(Tx_name)

# Error rates for all founder groups without missing founder
#######################################################################
# Calculate number of erroneous calls in founders according to 
# parsimonious call in decendent clones
# flsHom - at least 4 het calls in end-point clones
# flsHet - at least 6 Hom calls and no het calls in end-point clones

all_error_rates <- SNPs_merge_finalGT %>% 
  filter(!Line %in% noAncestor, 
         QUAL_BYcall >= 1000, QUAL_RMcall >= 1000, !Cut) %>%
  errorFromPhylo(flsHom_support = 4, flsHet_support = 6, output_POSi = F)

overall_F_Hom_rate <- sum(all_error_rates$n_F_Hom)/sum(all_error_rates$n_Het_q)
overall_F_Het_rate <- sum(all_error_rates$n_F_Het)/sum(all_error_rates$n_Hom_q)
overall_F_BY_rate <- sum(all_error_rates$nRef)/sum(all_error_rates$n_Het_q)
overall_F_RM_rate <- sum(all_error_rates$nAlt)/sum(all_error_rates$n_Het_q)

all_detected_errors_POSi <- SNPs_merge_finalGT %>% 
  filter(!Line %in% noAncestor, 
         QUAL_BYcall >= 1000, QUAL_RMcall >= 1000, !Cut) %>%
  errorFromPhylo(flsHom_support = 4, flsHet_support = 6, output_POSi = T) %>% 
  filter(is_error)

detected_errors_POSi_counts <- all_detected_errors_POSi %>% count(POSi)
detected_errors_POSi_counts %>% ggplot() + geom_histogram(aes(x = n), binwidth = 1)
all_detected_errors_POSi %>% filter(POSi < 203500 | POSi > 203700)

mixed_anc_GT_POSi <- sitewise_GTs %>% 
  filter(!(fHet_anc >= 0.875 | fHet_anc <= 0.125)) %>% pull(POSi)

repeat_error_POSi <- all_detected_errors_POSi %>% 
  count(POSi) %>% filter(n > 1) %>% pull(POSi)

high_error_POSi <- c(mixed_anc_GT_POSi, repeat_error_POSi)

# SNPs_merge_finalGT <- SNPs_merge_finalGT %>% 
#   filter(!POSi %in% high_error_POSi)

# Identify and count markers that are heterozygous across nearly all founders
n_markers <- nrow(sitewise_GTs)
sitewise_GTs %>% 
  filter(!POSi %in% high_error_POSi) %>% nrow()
clean_markers <- sitewise_GTs %>% 
  filter(nHom_anc <= 1, nHet_anc >= 4,
         !POSi %in% high_error_POSi) %>% pull(POSi)
n_markers_clean <- length(clean_markers)

# Get sites that are het in the founders for LOH detection
#######################################################################

final_cols <- c("CHROM", "POS", "POSi", "Tx_name", "Tx_ID", "Tx", "Line", "Rep", "ID", 
                "GT", "GQ", "Ref_DP_final", "Alt_DP_final", "Sum_DP_final",
                "existing_SNP", "Cut")

# For each founder group, collect calls in all clones that are Het in the founder
LOH_SNPs <- SNPs_merge_finalGT %>% 
  ungroup() %>%
  filter(!POSi %in% high_error_POSi, 
         QUAL_BYcall >= 1000, QUAL_RMcall >= 1000, 
         !Cut) %>% 
  select(all_of(final_cols)) %>% 
  anc_GT_fltr(anc_GT = "0/1")
LOH_SNPs$Line <- droplevels(LOH_SNPs$Line)
LOH_SNPs$ID <- droplevels(LOH_SNPs$ID)

# Make vectors of all sites for all clones and all LOH sites for all clones
LOH_rows <- paste0(LOH_SNPs$ID, "_", LOH_SNPs$POSi)
SNP_rows <- paste0(SNPs_merge_finalGT$ID, "_", SNPs_merge_finalGT$POSi)
loh_rows <- SNP_rows %in% LOH_rows

# Find sites for each clone not included in the LOH set. 
# After excluding annotated SNPs, the remainder are putative de novo mutations
denovo_rows <- !SNP_rows %in% LOH_rows
denovo_SNPs <- SNPs_merge_finalGT %>% 
  ungroup() %>%
  filter(denovo_rows, !existing_SNP, !POSi %in% high_error_POSi,
         QUAL_BYcall >= 100, QUAL_RMcall >= 100, !Cut) %>% 
  select(all_of(c(final_cols, "REF", "ALT"))) 

# Write table files
# save(SNPs_merge_finalGT, file = SNPs_merge_filename)
# save(LOH_SNPs, file = LOH_SNPs_file)
# save(denovo_SNPs, file = denovo_SNPs_file)

load(file = LOH_SNPs_file)

# Some sites may still be prone to erroneous calls and give false LOH events
# These are detected as exhibiting extreme allele bias in their conversion
# We filter out those markers that have been converted at least 10 times
# toward one allele out of at least 12 marker samplings

false_LOH_markers <- LOH_SNPs %>% 
  filter(GT != "0/1") %>%
  count(POSi, CHROM, POS, GT) %>% 
  pivot_wider(names_from = GT, values_from = n) %>% 
  rename(Ref = `0/0`, Alt = `1/1`) %>%
  replace(is.na(.), 0) %>% 
  filter((Ref < 2 & Alt > 10) | (Alt < 2 & Ref > 10)) %>% # data.frame()
  pull(POSi)

# LOH_SNPs %>% 
#   filter(GT != "0/1") %>% 
#   count(POSi, CHROM, POS, GT) %>% 
#   pivot_wider(names_from = GT, values_from = n) %>% 
#   rename(Ref = `0/0`, Alt = `1/1`) %>%
#   replace(is.na(.), 0) %>% 
#   filter(POSi %in% false_LOH_markers) %>%
#   mutate(f_Ref = Ref/(Ref+Alt)) %>%
#   filter(Ref+Alt > 10) %>% pull(Ref) %>% min()
# ggplot() + geom_histogram(aes(f_Ref), binwidth = 0.025)

LOH_SNPs <- LOH_SNPs %>% 
  filter(!POSi %in% false_LOH_markers)

# Calculate number of sites sampled for rate statistics
###############################################################################
n_SNP_sites <- SNPs_merge_finalGT %>% 
  filter(!POSi %in% high_error_POSi, !Cut) %>% 
  distinct(POSi) %>% nrow()
n_LOH_sites <- LOH_SNPs %>% distinct(POSi) %>% nrow()
n_DN_sites <- denovo_SNPs %>% distinct(POSi) %>% nrow()
n_dropped_sites <- n_sites_mean - n_SNP_sites
n_dropped_DN <- n_dropped_sites + n_LOH_sites


# Different clones have different sets of markers. In order to estimate marker
# distance, we take the mean for each marker accross clones. 
marker_dist <- LOH_SNPs %>% 
  group_by(ID) %>%
  mutate(d_POS = POSi - lag(POSi)) %>% 
  group_by(POSi) %>% 
  summarize(m_d_POSi = mean(d_POS))
  # slice_sample(n = 1) %>%
  # select(POSi, d_POS)

marker_dist %>%
  ungroup() %>%
  summarize(mn = mean(m_d_POSi, na.rm = T),
            md = median(m_d_POSi, na.rm = T))

marker_dist %>%
  # filter(d_POS < 10000) %>%
  # slice_sample(n = 1000) %>%
  ggplot() + geom_histogram(aes(x = log10(m_d_POSi)), binwidth = 1/3) +
  scale_y_log10()


all_chr_DP <- LOH_SNPs %>%
  dplyr::group_by(Tx, Line, Rep, ID, CHROM) %>% 
  dplyr::summarise(Ref_DP = mean(Ref_DP_final, na.rm = T),
                   Alt_DP = mean(Alt_DP_final, na.rm = T),
                   mean_DP = mean(Sum_DP_final, na.rm = T),
                   f_Het = sum(GT == "0/1")/sum(GT != "./."))

all_chr_DP <- all_chr_DP %>% 
  group_by(ID) %>% 
  mutate(f_Ref_DP = Ref_DP/mean(mean_DP),
         f_Alt_DP = Alt_DP/mean(mean_DP),
         f_Sum_DP = mean_DP/mean(mean_DP))


# Calculate the fractional difference between the fractional depth of each Chr in each clone
# and the mean depth of that chromosome
all_chr_DP <- all_chr_DP %>% group_by(CHROM) %>% 
  mutate(f_CHROM_Ref_DP = f_Ref_DP/mean(f_Ref_DP),
         f_CHROM_Alt_DP = f_Alt_DP/mean(f_Alt_DP),
         f_CHROM_Sum_DP = f_Sum_DP/mean(f_Sum_DP))


ID_CHROM_aneuploid <- all_chr_DP %>%
  filter(((f_CHROM_Ref_DP > 1.8 | f_CHROM_Alt_DP > 1.8) & f_Het > 0.9) | 
           ((f_CHROM_Ref_DP < 0.2 | f_CHROM_Alt_DP < 0.2) & f_Het < 0.1)) %>%
  # select(ID, f_Het, contains("DP"))
  mutate(ID_CHROM = paste0(ID, "_", CHROM))


# Tables of LOH boundaries and counts for each clone ------
all_GT_bounds <- LOH_SNPs %>% 
  filter(Rep != "00") %>% 
  filter(POSi %in% clean_markers) %>%
  EstDataBounds(., rm_noData = T)

all_GT_bounds_aneuCorr <- all_GT_bounds %>% 
  mutate(ID_CHROM = paste0(ID, "_", CHROM)) %>% 
  filter(!ID_CHROM %in% ID_CHROM_aneuploid$ID_CHROM)

# Mark LOH regions belonging to complex events
all_GT_bounds_aneuCorr <- all_GT_bounds_aneuCorr %>% 
  MarkLOHcomplexes(., gap = 10000)

all_GT_bounds_merge <- all_GT_bounds_aneuCorr %>%
  MergeComplexLOHs() %>% 
  MarkTerminalLOHs(ancHet = LOH_SNPs)
all_GT_bounds_merge$Tx_ID <- Recode_Tx_ID(all_GT_bounds_merge$Tx)

n_conv_markers <- LOH_SNPs %>% filter(GT != "0/1") %>% nrow()
overall_n_F_Hom <- round(overall_F_Hom_rate * 
  sum(all_GT_bounds_aneuCorr$length))
n_singleton_LOH <- all_GT_bounds_merge %>% filter(GT != "0/1", length == 1) %>% nrow()
n_doubleton_LOH <- all_GT_bounds_merge %>% filter(GT != "0/1", length == 2) %>% nrow()
n_F_Hom_rem <- overall_n_F_Hom - n_singleton_LOH

doubletons <- data.frame(ID = paste0("c_", 1:n_doubleton_LOH), n = 2)
n_reps <- 1000
tot_stats <- data.frame(NULL)
while(n_reps > 0) {
  ID_sample <- doubletons %>% slice_sample(n = n_F_Hom_rem, replace = T)
  ID_counts <- ID_sample %>% count(ID)
  if(max(ID_counts$n) < 3) {
    rep_stats_t <- ID_counts %>% count(n, name = "n_n") %>% arrange(n)
    rep_stats <- data.frame(rep = n_reps, single = rep_stats_t[1, 2], error = rep_stats_t[2, 2])
    tot_stats <- rbind(tot_stats, rep_stats)
    n_reps <- n_reps - 1
  }
}

# tot_stats_2 <- tot_stats
mean(tot_stats$single)
mean(tot_stats$error)

# No longer use this error correction. Now, remove all
# single-marker LOH events.
# all_GT_bounds_merge <- all_GT_bounds_merge %>% 
#   MarkErrorLOH_pooled(., error_rate = overall_F_Hom_rate)

all_LOHbounds_merge_NS <- all_GT_bounds_merge %>% 
  filter(GT != "0/1", length > 1)

all_LOHcounts_merge_NS <- all_LOHbounds_merge_NS %>% 
  CountLOHevents(omitError = F)
all_LOHcounts_merge_NS <- CategoriesFromID(all_LOHcounts_merge_NS)

# all_LOHcounts_merge_EC <- all_GT_bounds_merge %>% 
#   CountLOHevents(omitError = T)
# all_LOHcounts_merge_EC <- CategoriesFromID(all_LOHcounts_merge_EC)

all_LOHcounts_merge <- all_GT_bounds_merge %>% 
  CountLOHevents(omitError = F)
all_LOHcounts_merge <- CategoriesFromID(all_LOHcounts_merge)

# rm(SNPs_merge_finalGT)

# write LOH tables
###############################################################################
# save(all_LOHbounds_merge_NS, file = bounds_filename)
# save(all_LOHcounts_merge_NS, file = countsEC_filename)
# save(all_LOHcounts_merge, file = counts_filename)


