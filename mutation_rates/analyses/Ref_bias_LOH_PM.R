# Check the effect of reference sequence on allele bias in calling genotypes and detecting LOH events
# Start with SNPs_merge from 02_mainDataFrames

final_dubHet_rate <- SNPs_merge_finalGT %>% 
  filter(GT == "0/1") %>% 
  count(Cut) %>% 
  pivot_wider(names_from = "Cut", values_from = "n", names_prefix = "Cut_") %>%
  mutate(f_Cut = Cut_TRUE/(Cut_TRUE + Cut_FALSE))

load(file = SNPs_merge_raw_filename)
load(file = LOH_SNPs_file)
# Create analog of SNPs_merge_finalGT with only BY call data
SNPs_merge_BY <- SNPs_merge %>% 
  filter(!(ID %in% c(contaminated, bad_seq, whGnm_aneu_ID) | 
             Line %in% whGnm_aneu_line)) %>%
  filter(!is.na(GT_BYcall), QUAL_BYcall >= 1000, GQ_BYcall >= 50) %>%
  select(!(contains("_RMcall") |
         contains("REF", ignore.case = F) |
         contains("ALT", ignore.case = F))) %>%
  mutate(Sum_DP_final = Ref_DP_BYcall + Alt_DP_BYcall)

colnames(SNPs_merge_BY) <- gsub("_BYcall", "_final", colnames(SNPs_merge_BY))
colnames(SNPs_merge_BY) <- gsub("GT_final", "GT", colnames(SNPs_merge_BY))
colnames(SNPs_merge_BY) <- gsub("GQ_final", "GQ", colnames(SNPs_merge_BY))


# Fill in information, standardize column names, drop levels, 
# find annotated SNPs from RMxBY tables
SNPs_merge_BY <- SNPs_merge_BY %>% FormatSNPsMerge()

# Remove all calls with insufficient support
# Het - at least three reads for each allele
# Hom - at least three reads supporting allele
SNPs_merge_BY <- SNPs_merge_BY %>%
  filter(!(GT == "0/1" & (Ref_DP_final <= 3 | Alt_DP_final <= 3)),
         !(GT == "1/1" & Alt_DP_final <= 3),
         !(GT == "0/0" & Ref_DP_final <= 3)) 


# Mark heterozygous calls that have excessive allele imbalance
SNPs_merge_BY$Cut <- MarkDubiousHets(SNPs_merge_BY, .alpha = 0.05)

BYcall_dubHet_rate <- SNPs_merge_BY %>% 
  filter(GT == "0/1") %>% 
  count(Cut) %>% 
  pivot_wider(names_from = "Cut", values_from = "n", names_prefix = "Cut_") %>%
  mutate(f_Cut = Cut_TRUE/(Cut_TRUE + Cut_FALSE))

# Impute genotypes of missing ancestors from end-point clone data
SNPs_merge_BY <- SNPs_merge_BY %>% 
  InferMissingFounders(missing_founders = noAncestor)
SNPs_merge_BY <- SNPs_merge_BY %>% FormatSNPsMerge()

# SNPs_merge_BY %>% distinct(ID, .keep_all = T) %>% count(Line)
# save(SNPs_merge_BY, file = "./data/int/SNPs_merge_BYcall_2022_03.RData")

# Count each genotype for each site across founder clones and 
# end-point clones independently
sitewise_GTs_BYcall <- SNPs_merge_BY %>% 
  filter(!Cut) %>%
  site_genotype_stats(group = "all") %>%
  mutate(nHom_anc = nRef_anc + nRef_anc)

all_error_rates_BYcall <- SNPs_merge_BY %>% 
  filter(!Line %in% noAncestor, !Cut) %>%
  errorFromPhylo(flsHom_support = 4, flsHet_support = 6, output_POSi = F)

F_Hom_rate_BYcall <- sum(all_error_rates_BYcall$n_F_Hom)/sum(all_error_rates_BYcall$n_Het_q)
F_Het_rate_BYcall <- sum(all_error_rates_BYcall$n_F_Het)/sum(all_error_rates_BYcall$n_Hom_q)
F_BY_rate_BYcall <- sum(all_error_rates_BYcall$nRef)/sum(all_error_rates_BYcall$n_Het_q)
F_RM_rate_BYcall <- sum(all_error_rates_BYcall$nAlt)/sum(all_error_rates_BYcall$n_Het_q)

clean_markers <- sitewise_GTs_BYcall %>% 
  filter(nHom_anc <= 1, nHet_anc >= 4) %>% pull(POSi)

final_cols <- c("CHROM", "POS", "POSi", "Tx_name", "Tx_ID", "Tx", "Line", "Rep", "ID", 
                "GT", "GQ", "Ref_DP_final", "Alt_DP_final", "Sum_DP_final",
                "existing_SNP", "Cut")

# For each founder group, collect calls in all clones that are Het in the founder
LOH_SNPs_BYcall <- SNPs_merge_BY %>% 
  ungroup() %>%
  filter(!Cut) %>% 
  select(all_of(final_cols)) %>% 
  anc_GT_fltr(anc_GT = "0/1")


LOH_SNPs_BYcall$Line <- droplevels(LOH_SNPs_BYcall$Line)
LOH_SNPs_BYcall$ID <- droplevels(LOH_SNPs_BYcall$ID)

# Get GT tract boundaries
all_GT_bounds_BYcall <- LOH_SNPs_BYcall %>% 
  filter(Rep != "00") %>% 
  filter(POSi %in% clean_markers) %>%
  EstDataBounds(., rm_noData = T)

# Mark LOH regions belonging to complex events
all_GT_bounds_merge_BYcall <- all_GT_bounds_BYcall %>% 
  MarkLOHcomplexes(., gap = 10000) %>%
  MergeComplexLOHs() %>% 
  MarkErrorLOH(., error_rate = F_Hom_rate_BYcall) %>% 
  MarkTerminalLOHs(ancHet = LOH_SNPs_BYcall)

all_LOHcounts_merge_EC_BYcall <- all_GT_bounds_merge_BYcall %>% 
  CountLOHevents(omitError = T)
all_LOHcounts_merge_EC_BYcall <- CategoriesFromID(all_LOHcounts_merge_EC_BYcall)

all_LOHcounts_merge_BYcall <- all_GT_bounds_merge_BYcall %>% 
  CountLOHevents(omitError = F)
all_LOHcounts_merge_BYcall <- CategoriesFromID(all_LOHcounts_merge_BYcall)

all_LOHcounts_merge_BYcall <- all_LOHcounts_merge_BYcall %>% 
  mutate(n_BY = n_BYsmpl + n_BYcmplx,
         n_RM = n_RMsmpl + n_RMcmplx,
         callset = "BY")


# Create analog of SNPs_merge_finalGT with only RM call data
SNPs_merge_RM <- SNPs_merge %>% 
  filter(!(ID %in% c(contaminated, bad_seq, whGnm_aneu_ID) | 
             Line %in% whGnm_aneu_line)) %>%
  filter(!is.na(GT_RMcall), QUAL_RMcall >= 1000, GQ_RMcall >= 50) %>%
  select(!(contains("_BYcall") |
             contains("REF", ignore.case = F) |
             contains("ALT", ignore.case = F))) %>%
  mutate(Sum_DP_final = Ref_DP_RMcall + Alt_DP_RMcall)

colnames(SNPs_merge_RM) <- gsub("_RMcall", "_final", colnames(SNPs_merge_RM))
colnames(SNPs_merge_RM) <- gsub("GT_final", "GT", colnames(SNPs_merge_RM))
colnames(SNPs_merge_RM) <- gsub("GQ_final", "GQ", colnames(SNPs_merge_RM))


# Fill in information, standardize column names, drop levels, 
# find annotated SNPs from RMxBY tables
SNPs_merge_RM <- SNPs_merge_RM %>% FormatSNPsMerge()

# Remove all calls with insufficient support
# Het - at least three reads for each allele
# Hom - at least three reads supporting allele
SNPs_merge_RM <- SNPs_merge_RM %>%
  filter(!(GT == "0/1" & (Ref_DP_final <= 3 | Alt_DP_final <= 3)),
         !(GT == "1/1" & Alt_DP_final <= 3),
         !(GT == "0/0" & Ref_DP_final <= 3)) 


# Mark heterozygous calls that have excessive allele imbalance
SNPs_merge_RM$Cut <- MarkDubiousHets(SNPs_merge_RM, .alpha = 0.05)

# Impute genotypes of missing ancestors from end-point clone data
SNPs_merge_RM <- SNPs_merge_RM %>% 
  InferMissingFounders(missing_founders = noAncestor)
SNPs_merge_RM <- SNPs_merge_RM %>% FormatSNPsMerge()

# SNPs_merge_RM %>% distinct(ID, .keep_all = T) %>% count(Line)
save(SNPs_merge_RM, file = "./data/int/SNPs_merge_RMcall_2022_03.RData")

RMcall_dubHet_rate <- SNPs_merge_RM %>% 
  filter(GT == "0/1") %>% 
  count(Cut) %>% 
  pivot_wider(names_from = "Cut", values_from = "n", names_prefix = "Cut_") %>%
  mutate(f_Cut = Cut_TRUE/(Cut_TRUE + Cut_FALSE))

# Count each genotype for each site across founder clones and 
# end-point clones independently
sitewise_GTs_RMcall <- SNPs_merge_RM %>% 
  filter(!Cut) %>%
  site_genotype_stats(group = "all") %>%
  mutate(nHom_anc = nRef_anc + nRef_anc)

all_error_rates_RMcall <- SNPs_merge_RM %>% 
  filter(!Line %in% noAncestor, !Cut) %>%
  errorFromPhylo(flsHom_support = 4, flsHet_support = 6, output_POSi = F)


F_Hom_rate_RMcall <- sum(all_error_rates_RMcall$n_F_Hom)/sum(all_error_rates_RMcall$n_Het_q)
F_Het_rate_RMcall <- sum(all_error_rates_RMcall$n_F_Het)/sum(all_error_rates_RMcall$n_Hom_q)
F_BY_rate_RMcall <- sum(all_error_rates_RMcall$nRef)/sum(all_error_rates_RMcall$n_Het_q)
F_RM_rate_RMcall <- sum(all_error_rates_RMcall$nAlt)/sum(all_error_rates_RMcall$n_Het_q)

clean_markers <- sitewise_GTs_RMcall %>% 
  filter(nHom_anc <= 1, nHet_anc >= 4) %>% pull(POSi)

# For each founder group, collect calls in all clones that are Het in the founder
LOH_SNPs_RMcall <- SNPs_merge_RM %>% 
  ungroup() %>%
  filter(!Cut) %>% 
  select(all_of(final_cols)) %>% 
  anc_GT_fltr(anc_GT = "0/1")
LOH_SNPs_RMcall$Line <- droplevels(LOH_SNPs_RMcall$Line)
LOH_SNPs_RMcall$ID <- droplevels(LOH_SNPs_RMcall$ID)

# Get GT tract boundaries
all_GT_bounds_RMcall <- LOH_SNPs_RMcall %>% 
  filter(Rep != "00") %>% 
  filter(POSi %in% clean_markers) %>%
  EstDataBounds(., rm_noData = T)

# Mark LOH regions belonging to complex events
all_GT_bounds_merge_RMcall <- all_GT_bounds_RMcall %>% 
  MarkLOHcomplexes(., gap = 10000) %>%
  MergeComplexLOHs() %>% 
  MarkErrorLOH(., error_rate = F_Hom_rate_RMcall) %>% 
  MarkTerminalLOHs(ancHet = LOH_SNPs_RMcall)

all_LOHcounts_merge_EC_RMcall <- all_GT_bounds_merge_RMcall %>% 
  CountLOHevents(omitError = T)
all_LOHcounts_merge_EC_RMcall <- CategoriesFromID(all_LOHcounts_merge_EC_RMcall)

all_LOHcounts_merge_RMcall <- all_GT_bounds_merge_RMcall %>% 
  CountLOHevents(omitError = F)
all_LOHcounts_merge_RMcall <- CategoriesFromID(all_LOHcounts_merge_RMcall)

all_LOHcounts_merge_RMcall <- all_LOHcounts_merge_RMcall %>% 
  mutate(n_BY = n_BYsmpl + n_BYcmplx,
         n_RM = n_RMsmpl + n_RMcmplx,
         callset = "RM") 

all_LOHcounts_merge <- all_LOHcounts_merge %>% 
  mutate(n_BY = n_BYsmpl + n_BYcmplx,
         n_RM = n_RMsmpl + n_RMcmplx,
         callset = "Final") 

LOHsums <- all_LOHcounts_merge %>% 
  summarize(n_BY_tot = sum(n_BY), n_RM_tot = sum(n_RM),
            f_BY = sum(n_BY)/(sum(n_BY) + sum(n_RM)))

LOHsums_BYcall <- all_LOHcounts_merge_BYcall %>% 
  summarize(n_BY_tot = sum(n_BY), n_RM_tot = sum(n_RM),
            f_BY = sum(n_BY)/(sum(n_BY) + sum(n_RM)))

LOHsums_RMcall <- all_LOHcounts_merge_RMcall %>% 
  summarize(n_BY_tot = sum(n_BY), n_RM_tot = sum(n_RM),
            f_BY = sum(n_BY)/(sum(n_BY) + sum(n_RM)))

loh_cols <- c("Tx_name", "Line", "ID", "n_LOH", "n_BY", "n_RM", "callset")
by_cols <- c("Tx_name", "Line", "ID")
ns <- c("", "_BYcall", "_RMcall")

all_LOHcounts_byCall <- rbind(all_LOHcounts_merge[, loh_cols], 
                              all_LOHcounts_merge_BYcall[, loh_cols], 
                              all_LOHcounts_merge_RMcall[, loh_cols])
all_LOHcounts_byCall$callset <- factor(all_LOHcounts_byCall$callset, 
                                       levels = c("Final", "BY", "RM"))
all_LOHcounts_byCall <- all_LOHcounts_byCall %>%
  mutate(f_BY = n_BY/(n_BY + n_RM))

all_LOHcounts_byCall_mean <- all_LOHcounts_byCall %>%
  group_by(callset) %>%
  summarize(f_BY_mean = mean(f_BY), 
            CI_BY_lo =  mean(f_BY) - se(f_BY)*1.96, 
            CI_BY_up =  mean(f_BY) + se(f_BY)*1.96,
            bi_CI_lo = 0.5 - sqrt(0.25/n())*1.96,
            bi_CI_up = 0.5 + sqrt(0.25/n())*1.96)

n_c <- max(n_clones_xTx$n)
bi_var <- n_c * 0.25
c_bi <- rbinom(100, n_c, 0.5)
c_bi_var <- var(c_bi)

f_bi_var <- bi_var/n_c^2
f_c_bi <- rbinom(100, n_c, 0.5)/100
f_c_bi_var <- var(f_c_bi)
f_c_bi_var
# N * p * (1 - p)


all_LOHcounts_byCall %>% 
  ggplot() + 
  geom_vline(data = all_LOHcounts_byCall_mean, aes(xintercept = CI_BY_lo),
             color = "grey20", size = 0.15) +
  geom_vline(data = all_LOHcounts_byCall_mean, aes(xintercept = CI_BY_up),
             color = "grey20", size = 0.15) +
  geom_rect(data = all_LOHcounts_byCall_mean, 
            aes(xmin = CI_BY_lo, xmax = CI_BY_up, ymin = -Inf, ymax = Inf), 
             fill = "grey20", alpha = 0.2) +
  geom_vline(data = all_LOHcounts_byCall_mean, aes(xintercept = f_BY_mean), 
             color = "grey20", size = 1) +
  geom_histogram(aes(x = n_BY/(n_BY + n_RM)), fill = "grey20", binwidth = 0.05) +
  xlim(0, 1) +
  xlab("Fraction of LOH events to BY allele") +
  ylab("Number of clones") +
  facet_grid(callset~.)

DP_bias_final <- LOH_SNPs %>% 
  filter(GT == "0/1") %>% 
  # ungroup() %>% 
  group_by(ID) %>%
  summarize(f_BY = sum(Ref_DP_final)/sum(Ref_DP_final + Alt_DP_final)) %>%
  mutate(callset = "Final")


DP_bias_BYcall <- LOH_SNPs_BYcall %>% filter(GT == "0/1") %>% 
  group_by(ID) %>%
  summarize(f_BY = sum(Ref_DP_final)/sum(Ref_DP_final + Alt_DP_final)) %>%
  mutate(callset = "BY")

DP_bias_RMcall <- LOH_SNPs_RMcall %>% 
  filter(GT == "0/1") %>% 
  group_by(ID) %>%
  summarize(f_BY = sum(Ref_DP_final)/sum(Ref_DP_final + Alt_DP_final)) %>%
  mutate(callset = "RM")

DP_bias_all <- rbind(DP_bias_final, DP_bias_BYcall, DP_bias_RMcall)
DP_bias_all$callset <- factor(DP_bias_all$callset, levels = c("Final", "BY", "RM"))

DP_bias_all_mean <- DP_bias_all %>%
  group_by(callset) %>%
  summarize(f_BY_mean = mean(f_BY), 
            CI_BY_lo =  mean(f_BY) - se(f_BY)*1.96, 
            CI_BY_up =  mean(f_BY) + se(f_BY)*1.96,
            bi_CI_lo = 0.5 - sqrt(0.25/n())*1.96,
            bi_CI_up = 0.5 + sqrt(0.25/n())*1.96)
DP_bias_all_mean$callset <- factor(DP_bias_all_mean$callset, levels = c("Final", "BY", "RM"))

DP_bias_all %>% 
  ggplot() + 
  geom_vline( aes(xintercept = 0.5),
             color = "grey10", size = 0.25) +
  geom_rect(data = DP_bias_all_mean,
            aes(xmin = CI_BY_lo, xmax = CI_BY_up, ymin = -Inf, ymax = Inf),
            fill = "grey30", alpha = 0.25) +
  geom_vline(data = DP_bias_all_mean, aes(xintercept = CI_BY_lo),
             color = "grey20", size = 0.15) +
  geom_vline(data = DP_bias_all_mean, aes(xintercept = CI_BY_up),
             color = "grey20", size = 0.15) +
  geom_vline(data = DP_bias_all_mean, aes(xintercept = f_BY_mean),
             color = "black", size = 0.5) +
  geom_histogram(aes(x = f_BY), fill = "grey30", binwidth = 0.002) +
  xlim(0.47, 0.53) +
  ylim(0, 100) +
  xlab("Fraction of heterozygous call reads to BY allele") +
  ylab("Number of clones") +
  facet_grid(callset~.) +
  theme(panel.grid.minor.y = element_blank(),
        text = element_text(size = 18))

DP_bias_all %>% 
  ggplot(aes(x = callset, y = f_BY)) + 
  geom_boxplot(width = 0.3) +
  geom_line(aes(x = callset, y = f_BY, group = ID), 
            color = "grey50", alpha = 0.5)  +
  geom_point() +
  # ylim(0.48, 0.52) +
  theme(panel.grid.minor = element_blank(),
        text = element_text(size = 18))

# bias_merge <- merge(DP_bias_all, 
#                     all_LOHcounts_byCall[, c("ID", "callset", "f_BY")], 
#                     by = c("ID", "callset"), suffixes = c("_DP", "_LOH"))
# bias_merge %>% 
#   ggplot() + 
#   geom_jitter(aes(x = f_BY_DP, y = f_BY_LOH, color = callset), width = 0.01, height = 0.02)

# Final Figure ################################################################

RMvsBY_DP_rstest <- wilcox.test(f_BY ~ callset, data = DP_bias_all %>% filter(callset != "Final"))
FvsBY_DP_rstest <- wilcox.test(f_BY ~ callset, data = DP_bias_all %>% filter(callset != "RM"))
RMvsF_DP_rstest <- wilcox.test(f_BY ~ callset, data = DP_bias_all %>% filter(callset != "BY"))
rstest_DP_list <- list(RMvsBY = RMvsBY_DP_rstest, FvsBY = FvsBY_DP_rstest, RMvsF = RMvsF_DP_rstest)
rstest_DP_df <- data.frame(test_pair = names(rstest_DP_list), 
                           callset1 = c("RM", "Final", "RM"),
                           callset2 = c("BY", "BY", "Final"),
                           p_value = unlist(lapply(rstest_DP_list, function(x) x$p.value)))
rstest_DP_df <- rstest_DP_df %>% BHcorrection(., p_col = "p_value")
rstest_DP_df$p_value_rnd <- formatC(rstest_DP_df$p_value, digits = 2, format = "e")
rstest_DP_df <- cbind(rstest_DP_df, 
                      mean1 = DP_bias_all_mean$f_BY_mean[c(3, 1, 3)], 
                      mean2 = DP_bias_all_mean$f_BY_mean[c(2, 2, 1)])
rstest_DP_df$sig <- ifelse(rstest_DP_df$rejectNull_BH == 1, "*", "")



y_max_DP <- 110
DP_bias_plot <- DP_bias_all %>% 
  ggplot() + 
  geom_vline(aes(xintercept = 0.5),
             color = "grey10", size = 0.15) +
  geom_rect(xmin = 0.47, xmax = 0.53, ymin = y_max_DP + 1, ymax = Inf, 
            fill = "white", color = "white") +
  geom_bracket(data = rstest_DP_df, aes(xmin = mean1, xmax = mean2), 
               label = rstest_DP_df$sig,
               y.position = y_max_DP + 7, step.increase = 2,
               tip.length = 0.25, label.size = 8, vjust = 0.5) +
  geom_segment(data = DP_bias_all_mean, 
               aes(x = CI_BY_lo, xend = CI_BY_up, y = y_max_DP, yend = y_max_DP,
                   color = callset),
               size = 0.75) +
  geom_point(data = DP_bias_all_mean, 
             aes(x = f_BY_mean, y = y_max_DP, color = callset),
             size = 3.5) +
  geom_histogram(aes(x = f_BY, color = callset, fill = callset), 
                 position = position_dodge2(preserve = "single", padding = 0.0002), 
                 alpha = 1, binwidth = 0.002) +
  xlim(0.47, 0.53) +
  scale_y_continuous(breaks = seq(0, y_max_DP - 1, 20)) +
  # ylim(0, 100) +
  xlab("Fraction of heterozygous call reads to BY allele") +
  ylab("Number of clones") +
  scale_fill_manual(name = "Call set", values = c("grey45", "#FFC857", "#0D9FEE")) +
  scale_color_manual(values = c("grey45", "#FFC857", "#0D9FEE"), guide = NULL) +
  # facet_grid(callset~.) +
  theme(panel.grid.minor.y = element_blank(),
        text = element_text(size = 18),
        legend.position = c(0.95, 0.85))

DP_bias_plot

RMvsBY_LOH_rstest <- wilcox.test(f_BY ~ callset, data = all_LOHcounts_byCall %>% filter(callset != "Final"))
FvsBY_LOH_rstest <- wilcox.test(f_BY ~ callset, data = all_LOHcounts_byCall %>% filter(callset != "RM"))
RMvsF_LOH_rstest <- wilcox.test(f_BY ~ callset, data = all_LOHcounts_byCall %>% filter(callset != "BY"))
rstest_LOH_list <- list(RMvsBY = RMvsBY_LOH_rstest, FvsBY = FvsBY_LOH_rstest, RMvsF = RMvsF_LOH_rstest)
rstest_LOH_df <- data.frame(test_pair = names(rstest_DP_list), 
                            callset1 = c("RM", "Final", "RM"),
                            callset2 = c("BY", "BY", "Final"), 
                            p_value = unlist(lapply(rstest_LOH_list, function(x) x$p.value)))
rstest_LOH_df <- rstest_LOH_df %>% BHcorrection(., p_col = "p_value")
rstest_LOH_df <- cbind(rstest_LOH_df, 
                      mean1 = all_LOHcounts_byCall_mean$f_BY_mean[c(3, 1, 3)], 
                      mean2 = all_LOHcounts_byCall_mean$f_BY_mean[c(2, 2, 1)])
rstest_LOH_df$sig <- ifelse(rstest_LOH_df$rejectNull_BH == 1, "*", "")

y_max_LOH <- 42
LOH_bias_plot <- all_LOHcounts_byCall %>% 
  ggplot() + 
  geom_rect(xmin = 0, xmax = 1, ymin = y_max_LOH + 1, ymax = Inf, fill = "white", color = "white") +
  geom_bracket(data = rstest_LOH_df, aes(xmin = mean1, xmax = mean2),
               label = rstest_LOH_df$sig,
               y.position = y_max_LOH + 3, step.increase = 1.5,
               tip.length = 0.25, label.size = 8, vjust = 0.5) +
  geom_segment(data = all_LOHcounts_byCall_mean, 
               aes(x = CI_BY_lo, xend = CI_BY_up, 
                   y = y_max_LOH, yend = y_max_LOH,
                   color = callset),
               size = 0.75) +
  geom_point(data = all_LOHcounts_byCall_mean, 
               aes(x = f_BY_mean, y = y_max_LOH,
                   color = callset),
               size = 3.5) +
  geom_histogram(aes(x = f_BY, fill = callset), 
                 position = position_dodge2(preserve = "single"), 
                 # color = "grey20", 
                 binwidth = 0.05) +
  scale_fill_manual(values = c("grey45", "#FFC857", "#0D9FEE")) +
  scale_color_manual(values = c("grey45", "#FFC857", "#0D9FEE")) +
  scale_y_continuous(breaks = seq(0, y_max_LOH -1, 10)) +
  xlim(0, 1) +
  xlab("Fraction of LOH events to BY homolog") +
  ylab("Number of clones") +
  # facet_grid(callset~.) +
  theme(panel.grid.minor.y = element_blank(),
        text = element_text(size = 18),
        legend.position = "none")

LOH_bias_plot

bias_figure <- plot_grid(DP_bias_plot, LOH_bias_plot,
                    labels = c("A", "B"),
                    align = "v",
                    scale = 1,
                    label_size = 20,
                    hjust = 0,
                    ncol = 1, nrow = 2)

bias_figure

ggsave(file.path(outIntDir, "Allele_bias_combo_2022_03.png"), 
       plot = bias_figure,
       device = "png",
       width = 11, height = 8.5, 
       units = "in",
       dpi = 600)
