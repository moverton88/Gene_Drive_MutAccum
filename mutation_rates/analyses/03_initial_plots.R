
# For each final callset, calculate depth statistics by clone
###############################################################################
BYcall_allele_depth <- SNPs_merge_finalGT %>% 
  filter(GT_BYcall == "0/1", GQ_BYcall >= 50, !(GT_BYcall == "0/1" & Alt_DP_BYcall <= 3)) %>% 
  ungroup() %>% group_by(ID) %>% 
  summarize(call = "BY", 
            tot_BY_DP = sum(Ref_DP_BYcall), 
            tot_RM_DP = sum(Alt_DP_BYcall),
            fBY_DP = sum(Ref_DP_BYcall)/(sum(Ref_DP_BYcall) + sum(Alt_DP_BYcall)))

RMcall_allele_depth <- SNPs_merge_finalGT %>% 
  filter(GT_RMcall == "0/1", GQ_RMcall >= 50, !(GT_RMcall == "0/1" & Alt_DP_RMcall <= 3)) %>% 
  ungroup() %>% group_by(ID) %>% 
  summarize(call = "RM", 
            tot_BY_DP = sum(Ref_DP_RMcall), 
            tot_RM_DP = sum(Alt_DP_RMcall),
            fBY_DP = sum(Ref_DP_RMcall)/(sum(Ref_DP_RMcall) + sum(Alt_DP_RMcall)))

merge_allele_depth <- SNPs_merge_finalGT %>% 
  filter(GT == "0/1", GQ >= 50, !(GT == "0/1" & Alt_DP_final <= 3)) %>% 
  ungroup() %>% group_by(ID) %>% 
  summarize(call = "Merge", 
            tot_BY_DP = sum(Ref_DP_final), 
            tot_RM_DP = sum(Alt_DP_final),
            fBY_DP = sum(Ref_DP_final)/(sum(Sum_DP_final)))

# final depths after removing dubious Het calls
final_allele_depth <- SNPs_merge_finalGT %>% 
  filter(!Cut, GT == "0/1", GQ >= 50, !(GT == "0/1" & Alt_DP_final <= 3)) %>% 
  ungroup() %>% group_by(ID) %>% 
  summarize(call = "Final", 
            tot_BY_DP = sum(Ref_DP_final), 
            tot_RM_DP = sum(Alt_DP_final),
            fBY_DP = sum(Ref_DP_final)/(sum(Sum_DP_final)))

allele_depth <- rbind(BYcall_allele_depth, RMcall_allele_depth, merge_allele_depth, final_allele_depth)
allele_depth <- allele_depth %>% mutate(total_DP = (tot_BY_DP + tot_RM_DP))
allele_depth$call <- factor(allele_depth$call, levels = c("BY", "RM", "Merge", "Final"))

fBY_BYcVsRMc_perm <- perm_test(allele_depth, cat_var = "call", cat_names = c("BY", "RM"), 
                               response_var = "fBY_DP", test_stat = mean, n_perms = 10000)
fBY_BYcVsRMc_perm <- c(fBY_BYcVsRMc_perm, cat_1 = "BY", cat_2 = "RM")

fBY_BYcVsMerge_perm <- perm_test(allele_depth, cat_var = "call", cat_names = c("BY", "Merge"), 
                                 response_var = "fBY_DP", test_stat = mean, n_perms = 10000)
fBY_BYcVsMerge_perm <- c(fBY_BYcVsMerge_perm, cat_1 = "BY", cat_2 = "Merge")

fBY_RMcVsMerge_perm <- perm_test(allele_depth, cat_var = "call", cat_names = c("RM", "Merge"), 
                                 response_var = "fBY_DP", test_stat = mean, n_perms = 10000)
fBY_RMcVsMerge_perm <- c(fBY_RMcVsMerge_perm, cat_1 = "RM", cat_2 = "Merge")

fBY_MergeVsfinal_perm <- perm_test(allele_depth, cat_var = "call", cat_names = c("Merge", "Final"), 
                                   response_var = "fBY_DP", test_stat = mean, n_perms = 10000)
fBY_MergeVsfinal_perm <- c(fBY_MergeVsfinal_perm, cat_1 = "Merge", cat_2 = "Final")

fBY_perm_merge <- as.data.frame(bind_rows(fBY_BYcVsRMc_perm, fBY_BYcVsMerge_perm, 
                                      fBY_RMcVsMerge_perm, fBY_MergeVsfinal_perm))
fBY_perm_merge$cat_1 <- factor(fBY_perm_merge$cat_1, levels = c("BY", "RM", "Merge", "Final"))
fBY_perm_merge$cat_2 <- factor(fBY_perm_merge$cat_2, levels = c("BY", "RM", "Merge", "Final"))

fBY_perm_merge <- cbind(fBY_perm_merge, sig = ifelse(fBY_perm_merge$rejectNull, "*", " "))
# allele_depth %>% filter(call == "BY") %>% pull("fBY_DP") %>% t.test(., mu = 0.5)
# allele_depth %>% filter(call == "RM") %>% pull("fBY_DP") %>% t.test(., mu = 0.5)
# allele_depth %>% filter(call == "Merge") %>% pull("fBY_DP") %>% t.test(., mu = 0.5)
# allele_depth %>% filter(call == "Final") %>% pull("fBY_DP") %>% t.test(., mu = 0.5)

# Plot frequency of BY reference allele among callsets colored by total depth
###############################################################################
DP_mid <- log10(mean(allele_depth$total_DP))
y_range <- 0.04
allele_mapping_bias_plot <- allele_depth %>% 
  ggplot() + 
  geom_hline(aes(yintercept = 0.5), size = 0.2) +
  geom_boxplot(aes(x = call, y = fBY_DP), outlier.colour = NA, width = 0.5) +
  geom_jitter(aes(x = call, y = fBY_DP, color = log10(total_DP)), 
              height = 0, width = 0.2, alpha = 0.7, size = 2) +
  geom_bracket(data = fBY_perm_merge %>% filter(rejectNull), aes(xmin = cat_1, xmax = cat_2, label = sig), 
               y.position = 0.5 - 0.01 + y_range, step.increase = 0.075, vjust = 0.7, label.size = 5, fontface = "bold") +
  ylim(0.5 - y_range, 0.5 + y_range) +
  scale_color_gradient2(low = "brown3", mid = "grey50", high = "blue1", midpoint = DP_mid) +
  xlab("Callset") + ylab("Frequency of reference allele") +
  theme(text = element_text(size = 16))

allele_mapping_bias_plot

ggsave(file.path(outIntDir, "allele_mapping_bias_plot_2021_11.png"), 
       plot = allele_mapping_bias_plot,
       device = "png",
       width = 10, height = 10, 
       units = "in",
       dpi = 600)

################

allele_freq <- SNPs_merge_finalGT %>% filter(GT == "0/1") %>%
  # sample_n(100000) %>%
  summarize(fRef_BY = Ref_DP_BYcall/(Ref_DP_BYcall + Alt_DP_BYcall),
            fRef_RM = Ref_DP_RMcall/(Ref_DP_RMcall + Alt_DP_RMcall),
            fRef_final = Ref_DP_final/(Ref_DP_final + Alt_DP_final),
            het_cut = Cut) 

allele_freq %>%
  ggplot() + 
  geom_histogram(aes(x = fRef_BY), fill = "brown2", color = "brown2",
                 binwidth = 0.025, alpha = 0.5) +
  geom_histogram(aes(x = fRef_RM), fill = "dodgerblue2", color = "dodgerblue2",
                 binwidth = 0.025, alpha = 0.5) +
  geom_histogram(aes(x = fRef_final), fill = "white", color = "grey20",
                 binwidth = 0.025, alpha = 0.3)

  
allele_depth %>% 
  ggplot() + geom_abline(color = "blue4") +
  geom_point(aes(x = tot_BY_DP/1E6, y = tot_RM_DP/1E6, shape = call, group = ID), size = 2) + 
  geom_line(aes(x = tot_BY_DP/1E6, y = tot_RM_DP/1E6, group = ID), color = "grey30", alpha = 0.5) +
  # geom_label_repel(data = subset(allele_depth, tot_BY_DP/tot_RM_DP < 0.8 & call == "RM"), 
  #                  aes(x = tot_BY_DP, y = tot_RM_DP, label = ID), label.size = 0, force = 1.3) + 
  xlab("Reads supporting BY allele") + ylab("Reads supporting RM allele") +
  scale_shape_manual(values = c(15, 16, 17, 18), name = "Reference")

allele_depth %>% 
  ggplot() + geom_histogram(aes(x = fBY_DP)) + facet_wrap(~call, ncol = 1)

SNPs_merge_finalGT %>% 
  filter(GT == "0/1", !Cut) %>%
  ggplot() + geom_histogram(aes(x = Ref_DP_final/(Ref_DP_final + Alt_DP_final)),
                            breaks = seq(0.025, 0.975, 0.05)) + 
  facet_wrap(~Line)

# Compare distribution of BY allele frequencies for het sites and dubious het sites
SNPs_merge_finalGT %>% 
  filter(GT == "0/1") %>%
  ggplot() + geom_histogram(aes(x = 1 - f_Alt), breaks = seq(0.025, 0.975, 0.05)) +
  facet_grid(Cut~., scales = "free_y")

# Error rate distribution
all_error_rates %>% ggplot() + geom_histogram(aes(x = F_Hom_rate), binwidth = 0.00005)

# Determine whether high LOH regions are due to error-prone genotyping by 
# comparing proportion of clones with Het GT, anc vs evo, at these sites
# and also agreement between callsets

sitewise_GTs %>% 
  filter(nHet_anc >= 7) %>% 
  ggplot() + 
  # geom_histogram(aes(nNA_evo))
  geom_point(aes(x = POS, y = nNA_evo), size = 0.25) +
  facet_wrap(~CHROM, scales = "free_x")


sitewise_GTs %>% 
  # filter(nNA_anc <= 10, nNA_evo <= 125) %>%
  filter(nHet_anc >= 7, fHet_anc > 0.9) %>%
  # filter(POSi > 8240000, POSi < 8300000) %>%
  filter(CHROM == "09") %>%
  arrange(fHet_evo) %>%
  head(30)

# SNPs_merge_finalGT %>% 
LOH_SNPs %>%
  filter(CHROM == "01") %>% 
  # filter(Rep != "00") %>%
  # filter(Tx_name == "WT") %>%
  # filter(POS > 14000, POS < 27700) %>%
  # filter(POSi > 400000) %>%
  # count(POSi) %>% 
  ggplot() + 
  # geom_point(aes(x = POSi, y = n))
  # geom_vline(aes(xintercept = 435670), color = "blue3") +
  geom_point(aes(x = POS, y = ID, color = GT), size = 0.4)
  
SNPs_merge_finalGT %>% 
  filter(CHROM == "01") %>% 
  filter(Rep == "00") %>%
  # filter(Tx_name == "WT") %>%
  filter(POS > 15000, POS < 27700) %>% count(POS, GT) %>% pivot_wider(names_from = GT, values_from = n)
  # filter(POS < 50000) %>% 
  distinct(POSi) %>% arrange(POSi)
  
sitewise_GTs %>% 
  # filter(nNA_anc <= 10, nNA_evo <= 125) %>%
  filter(CHROM == "09") %>% 
  # arrange(fHet_evo)
  filter(POSi >= 5792574, POSi <= 5832461) %>% tail(30)
  ggplot() + geom_point(aes(x = fHet_anc, y = fHet_evo)) +
  geom_vline(aes(xintercept = 0.9))

sitewise_GTs %>% 
  # filter(nNA_anc <= 10, nNA_evo <= 125) %>%
  filter(CHROM == "09") %>%
  filter(POSi >= 5792574, POSi <= 5832461) %>% 
  select(CHROM:POSi, nNA_anc, fHet_anc, nRef_evo:nNA_evo) %>% tail(50)

rand_markers <- SNPs_merge_finalGT %>% 
  filter(Rep == "00") %>% dplyr::distinct(POSi) %>% 
  slice_sample(n = 400) %>% pull(POSi)

SNPs_merge_finalGT %>% 
  filter(Rep != "00") %>% 
  filter(POSi == 5828248) %>%
  filter(GT != "0/1") %>%
  # filter(POSi >= 5792574, POSi <= 5832461) %>%
  # filter(POSi %in% rand_markers) %>% 
  select(CHROM:POSi, ID, GT_BYcall, GQ_BYcall, GT_RMcall, GQ_RMcall, GT, GQ) %>% View()
  dplyr::count(GT, GT_BYcall, GT_RMcall)

sitewise_GTs %>% 
  filter(fHet_anc < 0.9) %>%
  # filter(nNA_anc <= 10, nNA_evo <= 125) %>%
  # filter(CHROM == "09") %>%
  # filter(POSi > 5800000) %>%
  ggplot() + geom_point(aes(x = POSi, y = fHet_anc)) +
  facet_wrap(~CHROM, scales = "free_x")


sitewise_GTs %>% 
  # filter(nNA_anc <= 10, nNA_evo <= 125) %>%
  ggplot() + geom_histogram(aes(x = fHet_anc))

sitewise_GTs %>% 
  filter(fHet_anc > 0.9) %>%
  ggplot() + geom_histogram(aes(x = nNA_anc), binwidth = 1)


sitewise_GTs %>% 
  # filter(nNA_anc <= 10) %>%
  filter(fHet_anc > 0.9) %>%
  dplyr::count(fHet_anc > 0.9)
  dplyr::count(nNA_anc <= 10)
  
poor_GT_sites <- sitewise_GTs %>% 
  filter(fHet_anc < 0.9) %>%
  pull(POSi)

# Relationship of number of markers included in an LOH event and its estimated length
# Shows that LOHs supported by single marker have similar estimated length 
# distributions to those with more marker support

all_LOHbounds_merge_EC %>%
  filter(GT != "0/1", !is_error) %>%
  # filter(!is.na(gap)) %>%
  ggplot() + geom_jitter(aes(x = log10(length), y = log10(est_length)), height = 0.05, width = 0.05)


# Plot number of sites in each clone that have valid final genotypes ------

cv_plot_in <- ID_siteCounts
cvrgPlot <- ggplot() + 
  geom_jitter(data = subset(cv_plot_in, Rep != "00"), 
              aes(x = Tx_name, y = n, group = ID), size = 2, height = 0, width = 0.3, alpha = 0.8) + 
  geom_jitter(data=subset(cv_plot_in, Rep == "00"), 
              aes(x = Tx_name, y=n, group=ID), color="orange1", size = 2,  height = 0, width = 0.3) + 
  # geom_text_repel(data=subset(cv_plot_in, n < minSites), 
  #                 aes(x=Tx_name, y=n, label=ID, group=ID), point.padding=0.2, force=0.5, segment.colour = NA) +
  geom_hline(aes(yintercept = minSites), linetype="dashed") +
  # geom_rect(aes(xmin = 0, xmax = 4, ymin = 0, ymax = minSites), 
  #           fill = "grey30", alpha = 0.1) +
  xlab("Drive Construct") + ylab("Number of valid sites")

cvrgPlot
