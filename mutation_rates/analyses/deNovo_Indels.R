###############################################################################
# Indels

# Import BY indel callset and filter out repeats
###############################################################################
indels_BY <- HC_multiVcfParse_allVar(BY_indel_file, all_Alts = T) %>% 
  filter(Sum_DP >= 6) %>% 
  filter(QUAL >= 100)

# Use BED of repeat annotations to filter indels
i_BY_rpt <- rep(F, nrow(indels_BY))
for(i in 1:nrow(repeats_bed)) {
  # i <- 1
  s_i <- indels_BY$POSi >= repeats_bed$Start_POSi[i]
  e_i <- indels_BY$POSi <= repeats_bed$End_POSi[i]
  # s_i <- sapply(indels_BY$POSi, function(x) x >= repeats_bed$Start_POSi[i], USE.NAMES = F)
  # e_i <- sapply(indels_BY$POSi, function(x) x <= repeats_bed$End_POSi[i], USE.NAMES = F)
  i_BY_rpt <- i_BY_rpt | (s_i & e_i)
}

indels_BY <- indels_BY[!i_BY_rpt, ]

alleles_BY_raw <- indels_BY %>% select(REF:ALT6)
repeats_BY <- rep(F, nrow(alleles_BY_raw))
for(a in 1:length(alleles_BY_raw)) {
  # a <- 1
  a_BY <- alleles_BY_raw[, a]
  repeats_allele <- sapply(a_BY, function(x) find_indel_repeats(x), USE.NAMES = FALSE)
  repeats_BY <- repeats_BY | repeats_allele
}

indels_BY_valid <- indels_BY[!repeats_BY, ]
rm(indels_BY)

# valid_new_indels_BY_POSi <- indels_BY_valid[!repeats_BY, ] %>% filter(Rep == "00", GT == "0/0") %>% distinct(POSi) %>% pull(POSi)
# indels_BY_valid %>% filter(POSi %in% valid_new_indels_BY_POSi, Rep != "00", GT != "0/0") %>% nrow()

# Import RM callset, fix position indicies, and remove repeats
###############################################################################
indels_RM <- HC_multiVcfParse_allVar(RM_indel_file, all_Alts = T) %>% 
  filter(Sum_DP >= 6) %>% 
  filter(QUAL >= 100)

# Liftover positions BY positions to RM callset
indels_RM_lift <- indels_RM %>% 
  RMxBY_liftover(.,  chain_df = BYtoRMchainDf, lift_from = "RM", lift_to = "BY")

rm(indels_RM)


i_RM_rpt <- rep(F, nrow(indels_RM_lift))
for(i in 1:nrow(repeats_bed)) {
  # i <- 1
  s_i <- indels_RM_lift$POSi >= repeats_bed$Start_POSi[i]
  e_i <- indels_RM_lift$POSi <= repeats_bed$End_POSi[i]
  i_RM_rpt <- i_RM_rpt | (s_i & e_i)
}

indels_RM_lift <- indels_RM_lift[!i_RM_rpt, ]

alleles_RM_raw <- indels_RM_lift %>% select(REF:ALT6)
repeats_RM <- rep(F, nrow(alleles_RM_raw))
for(a in 1:length(alleles_RM_raw)) {
  # a <- 1
  a_RM <- alleles_RM_raw[, a]
  repeats_allele <- sapply(a_RM, function(x) find_indel_repeats(x), USE.NAMES = FALSE)
  repeats_RM <- repeats_RM | repeats_allele
}

indels_RM_valid <- indels_RM_lift[!repeats_RM, ]
rm(indels_RM_lift)

# Merge indel callsets by ID and position
###############################################################################
merge_cols <- c("CHROM", "POS", "POSi", "QUAL", "REF", paste0("ALT", 1:6), "GT", "GQ",
                "Ref_DP", paste0("Alt", 1:6, "_DP"),
                "Tx", "Line", "Rep", "ID")

by_cols <- c("CHROM", "POS", "POSi", "Tx", "Line", "Rep", "ID")

indels_merge <- merge(indels_BY_valid[, merge_cols], 
                      indels_RM_valid[, merge_cols], 
                    by = by_cols, all = T, 
                    suffixes = c("_BYcall", "_RMcall"), 
                    sort = F) %>% arrange(ID, POSi)

# Remove positions with annotated RMxBY polymorphisms
indels_merge <- indels_merge %>% 
  filter(!POSi %in% RMxBY_comp_vcf$POSi) %>% 
  filter(!POSi %in% RMxBY_indels_POSi)

# Remove sites where on of the calls is missing
indels_merge_valid <- indels_merge %>%
  filter(!is.na(GT_BYcall), GT_BYcall != "./.",
         !is.na(GT_RMcall), GT_RMcall != "./.")

# Get allele columns for each unique position
indels_merge_alleles <- indels_merge_valid %>% 
  distinct(POSi, .keep_all = T) %>% 
  select(POSi, REF_BYcall:ALT6_BYcall, REF_RMcall:ALT6_RMcall)

# Find sites wherein the number of alternative alleles
# is the same for both callsets. Only want sites with consistent
# alleles b/t callsets
n_BY_alleles <- apply(indels_merge_alleles %>% 
                        select(REF_BYcall:ALT6_BYcall), MARGIN = 1, function(x) 7 - sum(is.na(x)))

n_RM_alleles <- apply(indels_merge_alleles %>% 
                        select(REF_RMcall:ALT6_RMcall), MARGIN = 1, function(x) 7 - sum(is.na(x)))

i_same_n_alleles <- n_BY_alleles == n_RM_alleles

# Get only allele sets with same number of alleles
indels_alleles_same_n <- indels_merge_alleles[i_same_n_alleles, ]
# Filter entire callset to these positions
indels_merge_same_n <- indels_merge_valid %>% filter(POSi %in% indels_alleles_same_n$POSi)

alleles_BY <- indels_alleles_same_n %>% select(REF_BYcall:ALT6_BYcall)
alleles_RM <- indels_alleles_same_n %>% select(REF_RMcall:ALT6_RMcall)

# For each position, create a map of the index in the BY alleles that 
# corresponds to each RM allele
# creates a dataframe that maps where in the BY alleles the RM allele matches
match_i_df <- data.frame(POSi = indels_alleles_same_n$POSi)
for(a in 1:length(alleles_RM)) {
  # a <- 1
  a_RM <- alleles_RM[, a]
  a_match_m <- sapply(alleles_BY, function(x) x == a_RM)
  a_match_m[is.na(a_match_m)] <- F
  a_match_i <- apply(a_match_m, MARGIN = 1, function(x) if(sum(x) != 0) {which(x)} else {0})
  match_i_df <- cbind(match_i_df, a_match_i)
  colnames(match_i_df)[a + 1] <- colnames(alleles_BY)[a]
}

# Split GT columns into individual allele designations (0/1 -> 0 1)
GT_BY_num <- colsplit(indels_merge_same_n$GT_BYcall, "/", names = paste0(c("A_1", "A_2"), "_BYcall"))
GT_RM_num <- colsplit(indels_merge_same_n$GT_RMcall, "/", names = paste0(c("A_1", "A_2"), "_RMcall"))
indels_merge_same_n <- cbind(indels_merge_same_n, 
                             lapply(GT_BY_num, function(x) as.numeric(x)), 
                             lapply(GT_RM_num, function(x) as.numeric(x)))

# Add "final" allele columns to fix RM genotypes non-distructively
indels_merge_same_n$A_1_final <- indels_merge_same_n$A_1_BYcall
indels_merge_same_n$A_2_final <- indels_merge_same_n$A_2_BYcall
# indels_merge_match <- indels_merge_same_n %>% filter(GT_BYcall == GT_RMcall)
# indels_merge_mis <- indels_merge_same_n %>% filter(GT_BYcall != GT_RMcall)

# Match RM GT alleles to BY callset
indels_merge_same_n <- match_GT_alleles(match_df = match_i_df, indels_df = indels_merge_same_n)

# Remove calls where there is not a complete match in alleles
# between the callsets (yields a -1 in one of the final GT alleles)
# This drops GTs if there are not perfect allele matches, and so #!!!!!!!!!!!!#
# likely loses valid final GTs and even whole sites #!!!!!!!!!!!!!!!!!!!!!!!!!#
indels_merge_same_n_final <- indels_merge_same_n %>% 
  filter(!(A_1_final == -1 | A_2_final == -1)) 

# Get final RM genotypes from matched GT alleles
indels_merge_final <- indels_merge_same_n_final %>% 
  mutate(GT_RMcall = paste0(A_1_final, "/", A_2_final)) %>% 
  select(CHROM:ID, GT_BYcall:GQ_BYcall, GT_RMcall, GQ_RMcall)

# Reconcile GTs based on GQ, as is done in the SNP analysis
indels_merge_finalGT <- GenotypeFromGQ(indels_merge_final, naThrsh = 100, 
                                       include_DP = F, check_alleles = F) %>%
  filter(finalGT != "./.") %>% mutate(GT = finalGT, GQ = finalGQ)

# For each position, count number of distinct GT in the founders
indel_anc_GTs <- indels_merge_finalGT %>% 
  filter(Rep == "00") %>% count(POSi, GT) %>% 
  pivot_wider(names_from = GT, values_from = n, names_prefix = "GT_")
indel_anc_GTs[is.na(indel_anc_GTs)] <- 0

# Find sites for which there are no founders with an existing indel
# and get the positions
i_anc_noIndel <- apply(indel_anc_GTs[, -c(1,2)], 
                       MARGIN = 1, function(x) sum(x, na.rm = T) == 0)
anc_noIndel_POSi <- indel_anc_GTs$POSi[i_anc_noIndel]

# For each position without an indel in the founders, count the number
# of distinct GT in the end-point clones
indel_evo_GTs <- indels_merge_finalGT %>% 
  filter(Rep != "00", POSi %in% anc_noIndel_POSi) %>% count(POSi, GT) %>% 
  pivot_wider(names_from = GT, values_from = n, names_prefix = "GT_")

# Find sites for which each non-Ref GT is represented by only
# one end-point clone. Allows for distinct indels to arise at one position
i_evo_Indel <- apply(indel_evo_GTs[, -c(1,2)], 
                     MARGIN = 1, function(x) sum(x > 1, na.rm = T) == 0)
evo_Indel_POSi <- indel_evo_GTs$POSi[i_evo_Indel]

# Filter indel set to valid sites as determined above
final_indel_set <- indels_merge_finalGT %>% 
  filter(POSi %in% evo_Indel_POSi, Rep != "00", GT != "0/0")

# Some indels are complex and represented multiple times, 
# merge indel events (filter out the redundants) within 100bp
# of each other
final_indel_set <- final_indel_set %>% 
  group_by(ID) %>%
  mutate(diff = POSi - lag(POSi),
         lead_diff = lead(POSi) - POSi)
final_indel_set$diff[is.na(final_indel_set$diff)] <- 12.2E7
final_indel_set$lead_diff[is.na(final_indel_set$lead_diff)] <- 12.2E7

final_indel_set_merge <- final_indel_set %>% filter(diff > 100 & lead_diff > 100)

# For indel counts, need to include clones with 0 indels
# Get IDs for all end-point clones from indels and SNPs
indel_IDs <- indels_merge %>% filter(Rep != "00") %>% distinct(ID) %>% pull(ID)
union_evo_IDs <- union(indel_IDs, evo_IDs)
final_indel_set_merge$ID <- factor(final_indel_set_merge$ID, levels = union_evo_IDs)

Indel_final_counts <- final_indel_set %>% count(ID) 
Indel_final_counts <- indels_merge_finalGT %>% 
  filter(POSi %in% anc_noIndel_POSi, Rep != "00", GT != "0/0") %>% count(ID) 

evo_zeros <- data.frame(ID = union_evo_IDs[!union_evo_IDs %in% Indel_final_counts$ID], n = 0)

Indel_final_counts <- rbind(evo_zeros, Indel_final_counts) %>% arrange(ID)

Indel_final_counts <- Indel_final_counts %>% CategoriesFromID()

Indel_final_counts %>% ggplot() + geom_histogram(aes(x = n), binwidth = 1) + facet_grid(Tx~.)
Indel_final_counts %>% group_by(Tx) %>% summarize(m = mean(n), s = sd(n))



###############################################################################

new_indels_merge <- indels_merge_valid %>% filter(!POSi %in% exist_indels_POSi)
# Need to calculate number of each GT at each position to find indels that are only
# in one or two clones
indel_POSi <- indels_merge_valid %>% distinct(POSi) %>% pull(POSi)

# For each position in the indel df, determines whether there is more than
# one GT represented by more than one clone. At real indels, expect
# all but one clone to be "0/1". This may be too conservative
i_unique_indels <- sapply(indel_POSi, function(x) indels_merge_valid %>% filter(POSi == x) %>% 
         count(GT_BYcall) %>% summarize(n_large = sum(n > 2) == 1) %>% pull(n_large))
unique_indels_POSi <- indel_POSi[i_unique_indels]

# Instead, we can get all founder genotypes and all end-point genotypes 
# at each site. Then filter to sites with at most one founder with a
# non 0/0 genotype and a max of 5 end-point clones with non 0/0 genotype
anc_BY_GT_wide <- indels_merge_valid %>% filter(Rep == "00", GT_BYcall != "./.") %>% distinct(POSi, ID, .keep_all = T) %>%
  select(POSi, ID, GT_BYcall) %>% pivot_wider(id_cols = c(POSi, ID), names_from = ID, values_from = GT_BYcall)

evo_BY_GT_wide <- indels_merge %>% filter(Rep != "00", GT_BYcall != "./.") %>% distinct(POSi, ID, .keep_all = T) %>%
  select(POSi, ID, GT_BYcall) %>% pivot_wider(id_cols = c(POSi, ID), names_from = ID, values_from = GT_BYcall)

anc_BY_GT_stats <- data.frame(POSi = anc_BY_GT_wide$POSi, 
                              n_hom_Ref = apply(anc_BY_GT_wide[, -1], MARGIN = 1, function(x) sum(x == "0/0", na.rm = T)),
                              n_non_Ref = apply(anc_BY_GT_wide[, -1], MARGIN = 1, function(x) sum(x != "0/0", na.rm = T)),
                              n_na = apply(anc_BY_GT_wide[, -1], MARGIN = 1, function(x) sum(is.na(x), na.rm = T)))

evo_BY_GT_stats <- data.frame(POSi = evo_BY_GT_wide$POSi, 
                              n_hom_Ref = apply(evo_BY_GT_wide[, -1], MARGIN = 1, function(x) sum(x == "0/0", na.rm = T)),
                              n_non_Ref = apply(evo_BY_GT_wide[, -1], MARGIN = 1, function(x) sum(x != "0/0", na.rm = T)),
                              n_na = apply(evo_BY_GT_wide[, -1], MARGIN = 1, function(x) sum(is.na(x), na.rm = T)))

anc_valid_POSi <- anc_BY_GT_stats %>% filter(n_non_Ref <= 1) %>% pull(POSi)

evo_valid_POSi <- evo_BY_GT_stats %>% 
  filter(n_non_Ref <= 5) %>% 
  # filter(n_hom_Ref >= 200, n_na <= 70) %>% 
  pull(POSi)

# Get the sites that are valid according to both the founder and e-p constraints
indel_valid_POSi <- intersect(anc_valid_POSi, evo_valid_POSi)

single_indel_POSi <- indels_merge_valid %>% 
  filter(POSi %in% indel_valid_POSi, Rep != "00", GT_BYcall != "0/0", GT_BYcall != "./.") %>% 
  count(POSi) %>% filter(n == 1) %>% pull(POSi)

multi_indel_POSi <- indels_merge_valid %>% 
  filter(POSi %in% indel_valid_POSi, Rep != "00", GT_BYcall != "0/0", GT_BYcall != "./.") %>% 
  count(POSi) %>% filter(n > 1) %>% pull(POSi)

indels_merge_valid %>% filter(POSi %in% single_indel_POSi, Rep != "00", GT_BYcall != "0/0", GT_BYcall != "./.") %>% count(Tx)
indels_merge_valid %>% filter(POSi %in% multi_indel_POSi, Rep != "00", GT_BYcall != "0/0", GT_BYcall != "./.") %>% count(Tx)


indels_merge_valid %>% # filter(POSi == 17289)
  # filter(POSi %in% anc_WT_POSi) %>% 
  filter(POSi %in% indel_valid_POSi) %>% 
  filter(Rep != "00", GT_BYcall != "0/0") %>% 
  select(POSi, ID, GT_BYcall) %>% arrange(POSi) %>% count(POSi, GT_BYcall)

evo_BY_GT_stats %>% filter(POSi %in% anc_valid_POSi)

indels_merge_valid %>% filter(POSi %in% unique_indels_POSi) %>% count(POSi)

indels_merge_valid %>% filter(POSi == 153968) %>% count(GT_BYcall)


# b_GT <- indels_merge_valid %>% select(c(GT_BYcall, GT_RMcall))

# GT_1 <- rep(".", nrow(alleles_RM_na_rm))
# GT_2 <- rep(".", nrow(alleles_RM_na_rm))

loose_match <- function(x, y) {
  if(is.na(x) | is.na(y)) {
    out <- F
  } else if(nchar(x) == nchar(y)) {
    out <- x == y
  } else if(nchar(x) > nchar(y)) {
    out <- length(grep(x, y)) > 0
  } else {
    out <- length(grep(y, x)) > 0
  }
  return(out)
}



new_single_indels_POSi <- new_indels_merge %>% 
  filter(Rep != "00", (GT_BYcall != "0/0" | GT_RMcall != "0/0")) %>% 
  count(POSi) %>% filter(n == 1) %>% pull(POSi)

new_indels_merge %>% 
  filter(Rep != "00", (GT_BYcall != "0/0" | GT_RMcall != "0/0"), 
         POSi %in% new_single_indels_POSi) %>% 
  count(Tx)


# Index putative biallelic SNPs by finding sites where the third allele is na
iBi_BY <- is.na(indels_BY$ALT2)
iBi_RM <- is.na(indels_RM_lift$ALT2)

# all_fltrd_RM_lift %>% filter(iBi_RM & iSNP_RM & !is.na(ALT2)) %>% distinct(POSi) %>% nrow

# Merge BY and RM dataframes by POSi and GT                        #
merge_cols <- c("CHROM", "POS", "POSi", "QUAL", "REF", "ALT1", "GT", "Ref_DP", "Alt1_DP",
                "GQ", "Tx", "Line", "Rep", "ID")

by_cols <- c("CHROM", "POS", "POSi", "Tx", "Line", "Rep", "ID")

indels_merge <- merge(indels_BY[iBi_BY, merge_cols], 
                      indels_RM[iBi_RM, merge_cols], 
                      by = by_cols, all = T, 
                      suffixes = c("_BYcall", "_RMcall"), 
                      sort = F) %>% arrange(ID, POSi)

colnames(indels_merge)[(grep("ALT1", colnames(indels_merge)))] <- c("ALT_BYcall", "ALT_RMcall")
colnames(indels_merge)[(grep("Alt1", colnames(indels_merge)))] <- c("Alt_DP_BYcall", "Alt_DP_RMcall")

indels_merge_swap <- SwapCrossedAlleles(SNP_df = indels_merge, annotated_POSi = RMxBY_indels_POSi, 
                                        to_match = "BYcall", to_swap = "RMcall")

indels_merge_swap <- indels_merge_swap %>% filter(!(is.na(bad_alleles)))

i_indels_match <- indels_merge_swap$REF_BYcall == indels_merge_swap$REF_RMcall & 
  indels_merge_swap$ALT_BYcall == indels_merge_swap$ALT_RMcall

indels_merge_swap <- indels_merge_swap[i_indels_match,]

indels_final <- indels_merge_swap %>% GenotypeFromGQ_mut() %>% filter(finalGT != "./.")

indels_final %>% filter(Rep == "00") %>% count(finalGT)
anc_homRef_POSi <- indels_final %>% filter(Rep == "00", finalGT == "0/0", finalGQ >= 50) %>% distinct(POSi) %>% pull(POSi) 

evo_indels <- indels_final %>% 
  filter(Rep != "00", POSi %in% anc_homRef_POSi) %>% 
  select(CHROM, POS, POSi, ID, finalGT) %>% pivot_wider(names_from = ID, values_from = finalGT)

evo_indels_n <- evo_indels %>% select(c("CHROM", "POS", "POSi"))

evo_indels_n$n_Het <- evo_indels %>% select(!c("CHROM", "POS", "POSi")) %>% apply(., MARGIN = 1, function(x) sum(x == "0/1", na.rm = T))
evo_indels_n$n_homAlt <- evo_indels %>% select(!c("CHROM", "POS", "POSi")) %>% apply(., MARGIN = 1, function(x) sum(x == "1/1", na.rm = T))
evo_indels_n$n_other <- evo_indels %>% 
  select(!c("CHROM", "POS", "POSi")) %>% 
  apply(., MARGIN = 1, function(x) sum(!x %in% c("0/0", "0/1", "1/1"), na.rm = T))


indels_final %>% count(finalGT)

evo_indels_n %>% ggplot() + geom_histogram(aes(x = n_homAlt))

valid_indels_POSi <- evo_indels %>% filter(n_Het == 1 & n_homAlt == 0 | n_homAlt == 1 & n_Het == 0) %>% distinct(POSi) %>% pull(POSi)

final_valid_indels <- indels_final %>% filter(POSi %in% valid_indels_POSi)
final_valid_indels %>% filter(Rep != "00", !finalGT == "0/0") %>% count(POSi)
final_valid_indels %>% filter(Rep != "00", !finalGT == "0/0") %>% arrange(POSi)


# de novo Indel analysis

for(l in levels(all_finalGT$Line)) {
  # l = "N_A"
  line_df <- vcf_df %>% filter(Line == l)
  anc_hom <- line_df %>% filter(Rep == "00" & GT == "0/0") %>% pull(POSi)
  line_df_wide <- line_df %>% 
    # filter(Rep != "00") %>% 
    select(POSi, ID, GT) %>% pivot_wider(names_from = ID, values_from = GT)
  line_sum_GTs <- line_df_wide %>% select(POSi, contains("00"))
  colnames(line_sum_GTs)[2] <- "ancGT"
  line_sum_GTs$n_Valid <- apply(line_df_wide[, -1], 1, function(x) sum(!is.na(x)))
  line_sum_GTs$n_refHom <- apply(line_df_wide[, -1], 1, function(x) sum(x == "0/0", na.rm = T))
  line_sum_GTs$n_Het <- apply(line_df_wide[, -1], 1, function(x) sum(x == "0/1", na.rm = T))
  line_sum_GTs$n_altHet <- apply(line_df_wide[, -1], 1, function(x) sum(x == "1/1", na.rm = T))
  valid_mut_POSi <- line_sum_GTs %>% 
    filter(n_Het == 1 & (ancGT == "0/0" | (is.na(ancGT) & n_refHom >= 2))) %>% 
    pull(POSi)
  mut_df_wide <- line_df_wide %>% filter(POSi %in% valid_mut_POSi)
  n_mut_ID <- apply(mut_df_wide[, -c(1)], 2, function(x) sum(x == "0/1" | x == "1/1", na.rm = T))
}

# Mark which loci contain indels
markIndels <- function(vcf_df) {
  # vcf_df <- all_fltrd_RM_POSadj
  alleleNames <- c("REF", colnames(vcf_df)[grep("ALT", colnames(vcf_df))])
  vcf_alleles <- vcf_df[, alleleNames]
  # nchar(all_fltrd_RM_POSadj[3,4])
  alleleLengths <- sapply(vcf_alleles, function(x) nchar(x)) 
  iIndels <- apply(alleleLengths, MARGIN = 1, function(x) max(x, na.rm = T)) > 1
  # vcf_df$Indel <- iIndels
  return(iIndels)
}

all_fltrd_RM_POSadj$Indel <- markIndels(all_fltrd_RM_POSadj)

# Counts the number of alleles at each locus that contain simple repeats
tallyRepeatAlleles <- function(vcf_df, n_HomoPs = 4, n_repeats = 4) {
  # vcf_df <- all_fltrd_RM_POSadj %>% filter(Indel == T) %>% slice(400:500)
  alleleNames <- c("REF", colnames(vcf_df)[grep("ALT", colnames(vcf_df))])
  # loops through allele columns and logs 0 for no homopolymer and 1 for 
  # presence of homopolymer. For each site, number of homopolymer alleles 
  # is summed
  sumHomoPs <- 0
  for(a in alleleNames) {
    # a <- "ALT1"
    iHomoPs <- lapply(vcf_df[, a], 
                      function(x) 
                        sum(rle(unlist(strsplit(x, split = "")))$lengths >= n_HomoPs))
    iHomoPs <- unlist(iHomoPs)
    is.na(iHomoPs) <- 0
    sumHomoPs <- sumHomoPs + iHomoPs
  }
  # sumAlleleRepeats <- 0
  sumRepeats <- 0
  for(a in alleleNames) {
    # Breaks down each allele string into chunks, counts runs of repeating patterns
    # k <- vcf_df[65, a]
    # allele_len <- nchar(k)
    # r = 3
    for(r in 2:3) {
      iRepeats <- lapply(vcf_df[, a], 
                         function(x) ifelse(nchar(x) > r*2,
                                            sum(rle(substring(x, 
                                                              seq.int(1, (nchar(x) - 1), r), 
                                                              seq.int(r, nchar(x), r))
                                            )$lengths >= n_repeats), 0))
      iRepeats <- unlist(iRepeats)
      iRepeats[is.na(iRepeats)] <- 0
      iRepeatsShft <- lapply(vcf_df[, a], 
                             function(x) ifelse(nchar(x) > r*2,
                                                sum(rle(substring(x, 
                                                                  seq.int(2, nchar(x), r), 
                                                                  seq.int(r + 1, nchar(x) + 1, r))
                                                )$lengths >= n_repeats), 
                                                0))
      iRepeatsShft <- unlist(iRepeatsShft)
      iRepeatsShft[is.na(iRepeatsShft)] <- 0
      if(r == 3) {
        iRepeatsShft2 <- lapply(vcf_df[, a], 
                                function(x) ifelse(nchar(x) > r*2,
                                                   sum(rle(substring(x, 
                                                                     seq.int(3, (nchar(x) + 1), r), 
                                                                     seq.int(r + 2, nchar(x) + 2, r))
                                                   )$lengths >= n_repeats), 0))
        iRepeatsShft2 <- unlist(iRepeatsShft2)
        iRepeatsShft2[is.na(iRepeatsShft2)] <- 0
      }
      sumRepeats <- sumRepeats + iRepeats + iRepeatsShft
    }
  }
  totalRepeats <- sumHomoPs + sumRepeats
  return(totalRepeats)
}


all_fltrd_RM_POSadj$n_repeats <- tallyRepeatAlleles(all_fltrd_RM_POSadj)

all_fltrd_RM_POSadj %>% filter(Indel == T, n_repeats == 0) %>% count(Tx)


markManyAlleleDepths <- function(vcf_df) {
  alleleDepthNames <- c("Ref_DP", colnames(vcf_df)[grep("Alt", colnames(vcf_df))])
  iMAD <- rowSums(vcf_df[, alleleDepthNames] != 0, na.rm = T)
}


all_fltrd_RM_POSadj$Rpts <- tallyRepeatAlleles(all_fltrd_RM_POSadj)
all_fltrd_RM_POSadj$MADs <- markManyAlleleDepths(all_fltrd_RM_POSadj)
