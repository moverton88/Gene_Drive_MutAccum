# Parallel process to 02_mainDataFrames.R,except with triploid founder groups
# F_B and F_E

triploid_SNPs_merge <- SNPs_merge %>% filter(Line %in% whGnm_aneu_line)
# save(triploid_SNPs_merge, file = triploid_merge_raw_filename)
# load(file = triploid_merge_raw_filename)

triploid_merge_finalGT <- triploid_SNPs_merge %>% 
  filter(!(ID %in% c(contaminated, bad_seq))) %>%
  GenotypeFromGQ(., baseThrsh = 50, naThrsh = 100, diffThrsh = 30) %>%
  filter(finalGT != "./.", finalGQ >= 50)

rm(triploid_SNPs_merge)

# Fill in information, standardize column names, drop levels, 
# find annotated SNPs from RMxBY tables
triploid_merge_finalGT <- triploid_merge_finalGT %>% FormatSNPsMerge()

# Remove all calls with insufficient support
# Het - at least three reads for each allele
# Hom - at least three reads supporting allele
triploid_merge_finalGT <- triploid_merge_finalGT %>%
  filter(!(GT == "0/1" & (Ref_DP_final <= 3 | Alt_DP_final <= 3)),
         !(GT == "1/1" & Alt_DP_final <= 3),
         !(GT == "0/0" & Ref_DP_final <= 3)) 

# Mark heterozygous calls that have excessive allele imbalance
triploid_merge_finalGT$Cut <- MarkDubiousHets(triploid_merge_finalGT, .alpha = 0.05, known_hets = clean_markers)


# Count each genotype for each site across founder clones and 
# end-point clones independently
triploid_sitewise_GTs <- triploid_merge_finalGT %>% 
  filter(!Cut) %>%
  site_genotype_stats(group = "all") %>%
  mutate(nHom_anc = nRef_anc + nRef_anc)

# Get IDs of all end-point clones with SNP data
triploid_all_IDs <- triploid_merge_finalGT %>% 
  distinct(ID, .keep_all = T) %>% 
  pull(ID) %>% as.character()

triploid_evo_IDs <- triploid_merge_finalGT %>% 
  distinct(ID, .keep_all = T) %>% 
  filter(Rep != "00") %>% pull(ID) %>% as.character()

triploid_n_founders <- triploid_merge_finalGT %>% 
  filter(Rep == "00") %>%
  distinct(ID, .keep_all = T) %>% count(Tx_name)

triploid_n_clones <- triploid_merge_finalGT %>% 
  filter(Rep != "00") %>% 
  pull(ID) %>% unique() %>% length()

triploid_n_clones_xTx <- triploid_merge_finalGT %>% 
  filter(Rep != "00") %>%
  distinct(ID, .keep_all = T) %>% count(Tx_name)

triploid_n_all_xTx <- triploid_merge_finalGT %>% 
  distinct(ID, .keep_all = T) %>% count(Tx_name)

# Error rates for all founder groups without missing founder
#######################################################################
# Calculate number of erroneous calls in founders according to 
# parsimonious call in decendent clones
# flsHom - at least 4 het calls in end-point clones
# flsHet - at least 6 Hom calls and no het calls in end-point clones

triploid_all_error_rates <- triploid_merge_finalGT %>% 
  filter(!Line %in% noAncestor, 
         QUAL_BYcall >= 1000, QUAL_RMcall >= 1000, !Cut) %>%
  errorFromPhylo(flsHom_support = 2, flsHet_support = 4, output_POSi = F)

# triploid_merge_finalGT <- triploid_merge_finalGT %>% 
#   filter(!POSi %in% high_error_POSi)

# Identify and count markers that are heterozygous across nearly all founders
triploid_n_markers <- sitewise_GTs %>% 
  filter(!POSi %in% high_error_POSi) %>% nrow()
# triploid_clean_markers <- triploid_sitewise_GTs %>% 
#   filter(nHom_anc == 0, nHet_anc >= 1,
#          !POSi %in% high_error_POSi) %>% pull(POSi)
# triploid_n_markers_clean <- length(triploid_clean_markers)

# Get sites that are het in the founders for LOH detection
#######################################################################

final_cols <- c("CHROM", "POS", "POSi", "Tx_name", "Tx_ID", "Tx", "Line", "Rep", "ID", 
                "GT", "GQ", "Ref_DP_final", "Alt_DP_final", "Sum_DP_final",
                "existing_SNP", "Cut")

# For each founder group, collect calls in all clones that are Het in the founder
triploid_LOH_SNPs <- triploid_merge_finalGT %>% 
  ungroup() %>%
  filter(QUAL_BYcall >= 1000, QUAL_RMcall >= 1000, 
         !Cut, !POSi %in% high_error_POSi) %>% 
  select(all_of(final_cols)) %>% 
  anc_GT_fltr(anc_GT = "0/1")
triploid_LOH_SNPs$Line <- droplevels(triploid_LOH_SNPs$Line)
triploid_LOH_SNPs$ID <- droplevels(triploid_LOH_SNPs$ID)

# Make vectors of all sites for all clones and all LOH sites for all clones
triploid_LOH_rows <- paste0(triploid_LOH_SNPs$ID, "_", triploid_LOH_SNPs$POSi)
triploid_SNP_rows <- paste0(triploid_merge_finalGT$ID, "_", triploid_merge_finalGT$POSi)
triploid_loh_rows <- triploid_SNP_rows %in% triploid_LOH_rows

# Find sites for each clone not included in the LOH set. 
# After excluding annotated SNPs, the remainder are putative de novo mutations
triploid_denovo_rows <- !triploid_SNP_rows %in% triploid_LOH_rows
triploid_denovo_SNPs <- triploid_merge_finalGT %>% 
  ungroup() %>%
  filter(triploid_denovo_rows, !existing_SNP,
         QUAL_BYcall >= 100, QUAL_RMcall >= 100, !Cut, 
         !POSi %in% high_error_POSi) %>% 
  select(all_of(final_cols)) 

# Write table files
# save(triploid_merge_finalGT, file = triploid_merge_filename)
# save(triploid_LOH_SNPs, file = triploid_LOH_SNPs_file)
# save(triploid_denovo_SNPs, file = triploid_denovo_SNPs_file)

triploid_dn_GT_vals <- triploid_denovo_SNPs %>% 
  site_genotype_stats()

# Get sitewise genotype stats for each founder group
triploid_dn_Line_GT_vals <- data.frame(NULL)
for(l in levels(triploid_denovo_SNPs$Line)) {
  # l <- "N_C"
  line_SNPs <- triploid_denovo_SNPs %>% filter(Line == l)
  line_SNPs$Line <- droplevels(line_SNPs$Line)
  line_GT_vals <- line_SNPs %>% site_genotype_stats()
  line_GT_vals$Line <- l
  triploid_dn_Line_GT_vals <- rbind(triploid_dn_Line_GT_vals, line_GT_vals)
}

# variant_other_lines <- dn_Line_GT_vals %>% 
#   filter(!(nHet_evo == 0 & nAlt_evo == 0)) %>% pull(POSi)
# Get sites for which there is a single polymorphism among all clones of a founder group
evo_single_mut_Line <- triploid_dn_Line_GT_vals %>% 
  filter(nRef_evo >= 4 & ((nHet_evo == 1 & nAlt_evo == 0) | (nAlt_evo == 1 & nHet_evo == 0))) %>% 
  select(Line, POSi, nRef_anc, nRef_evo, nHet_evo)

# Get positions were the above is true
evo_single_Line_POSi <- evo_single_mut_Line %>%
  distinct(POSi) %>% pull(POSi)

# Find sites for which there is a single polymorphism among all end-point clones and get positions
evo_single_mut_POSi <- triploid_dn_GT_vals %>% 
  filter(POSi %in% evo_single_Line_POSi & 
           !(POSi %in% evo_single_Line_POSi) &
           ((nHet_evo == 1 & nAlt_evo == 0) | 
              (nHet_evo == 0 & nAlt_evo == 1))) %>%
  pull(POSi)

# Sites for which all founders are homozygous reference
anc_hom_POSi <- triploid_dn_GT_vals %>% filter(fRef_anc == 1) %>% pull(POSi) %>% unique()

# evo_hom_POSi <- dn_GT_vals %>% filter(nRef_evo >= 30) %>% pull(POSi) %>% unique()

# Sites for which only a single end-point clone is heterozygous
evo_single_het_POSi <- triploid_dn_GT_vals %>% 
  filter(nHet_evo == 1 & nAlt_evo == 0) %>% pull(POSi) %>% unique()

# Sites for which only a single end-point clone is homozygous alt
evo_single_alt_POSi <- triploid_dn_GT_vals %>% 
  filter(nHet_evo == 0 & nAlt_evo == 1) %>% pull(POSi) %>% unique()

# Sites for which two end-point clones are heterozygous
evo_double_het_POSi <- triploid_dn_GT_vals %>% 
  filter(nHet_evo == 2 & nAlt_evo == 0) %>% pull(POSi) %>% unique()

# Unique heterozygous variant sites
anc_hom_mut_het_POSi <- intersect(anc_hom_POSi, evo_single_het_POSi)

# Unique homozygous variant sites
anc_hom_mut_alt_POSi <- intersect(anc_hom_POSi, evo_single_alt_POSi)

# Unique variant sites
triploid_mut_all_POSi <- union(union(anc_hom_mut_het_POSi, anc_hom_mut_alt_POSi), evo_single_mut_POSi)

triploid_mut_set <- triploid_denovo_SNPs %>% 
  filter(POSi %in% triploid_mut_all_POSi, Rep != "00", 
         GT %in% c("0/1", "1/1")) %>% arrange(ID, POSi)

triploid_mut_set <- triploid_mut_set %>% 
  group_by(ID) %>%
  mutate(diff = POSi - lag(POSi),
         lead_diff = lead(POSi) - POSi)
triploid_mut_set$diff[is.na(triploid_mut_set$diff)] <- 12.2E7
triploid_mut_set$lead_diff[is.na(triploid_mut_set$lead_diff)] <- 12.2E7

# triploid_mut_set %>% filter(diff < 100 | lead_diff < 100) %>% select(!QUAL_BYcall:Tx_name) %>% View()

# All clustered SNPs are in Ty LTRs and other types of repeats
triploid_mut_set_merge <- triploid_mut_set %>% filter(diff > 100 & lead_diff > 100)
triploid_mut_set_merge$CHROM <- factor(triploid_mut_set_merge$CHROM, levels = str_pad(1:16, width = 2, pad = "0"))
n_complex <- nrow(triploid_mut_set) - nrow(triploid_mut_set_merge)

triploid_n_SNMs <- nrow(triploid_mut_set_merge)

# Number of diploid positions
triploid_SNM_rate <- triploid_n_SNMs / (length(triploid_evo_IDs)*n_gens*mean_cover*2)
triploid_SNM_rate

###############################################################################
# Compare distribution of SNMs among clones
triploid_final_counts <- triploid_mut_set_merge %>% count(ID) 

triploid_evo_zeros <- data.frame(ID = evo_IDs[!evo_IDs %in% triploid_final_counts$ID], n = 0)

triploid_final_counts <- rbind(evo_zeros, triploid_final_counts) %>% arrange(ID)

triploid_final_counts <- triploid_final_counts %>% CategoriesFromID()


# triploid_final_counts %>% group_by(Line %in% noAncestor) %>% count(n)
triploid_mean_SNMrate <- triploid_final_counts %>% 
  ungroup() %>% group_by(Tx_name) %>% 
  summarise(total_SNMs = sum(n), mean_SNM = mean(n), sd_SNM = sd(n), n_clones = n())

triploid_mean_SNMrate$se_SNM <- triploid_mean_SNMrate$sd_SNM/sqrt(triploid_mean_SNMrate$n_clones)
triploid_mean_SNMrate$mean_rate <- triploid_mean_SNMrate$mean_SNM / n_gens / (mean_cover*2)
triploid_mean_SNMrate$se_rate <- triploid_mean_SNMrate$se_SNM / n_gens / (mean_cover*2)



