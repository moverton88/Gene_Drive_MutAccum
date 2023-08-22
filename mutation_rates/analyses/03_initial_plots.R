
# For each final callset, calculate depth statistics by clone
###############################################################################
P1call_allele_depth <- SNPs_merge_finalGT %>% 
  filter(GT_P1call == "0/1", GQ_P1call >= 50, !(GT_P1call == "0/1" & Alt_DP_P1call <= 3)) %>% 
  ungroup() %>% group_by(ID) %>% 
  summarize(call = "P1", 
            tot_P1_DP = sum(Ref_DP_P1call), 
            tot_P2_DP = sum(Alt_DP_P1call),
            fP1_DP = sum(Ref_DP_P1call)/(sum(Ref_DP_P1call) + sum(Alt_DP_P1call)))

P2call_allele_depth <- SNPs_merge_finalGT %>% 
  filter(GT_P2call == "0/1", GQ_P2call >= 50, !(GT_P2call == "0/1" & Alt_DP_P2call <= 3)) %>% 
  ungroup() %>% group_by(ID) %>% 
  summarize(call = "P2", 
            tot_P1_DP = sum(Ref_DP_P2call), 
            tot_P2_DP = sum(Alt_DP_P2call),
            fP1_DP = sum(Ref_DP_P2call)/(sum(Ref_DP_P2call) + sum(Alt_DP_P2call)))

merge_allele_depth <- SNPs_merge_finalGT %>% 
  filter(GT == "0/1", GQ >= 50, !(GT == "0/1" & Alt_DP_final <= 3)) %>% 
  ungroup() %>% group_by(ID) %>% 
  summarize(call = "Merge", 
            tot_P1_DP = sum(Ref_DP_final), 
            tot_P2_DP = sum(Alt_DP_final),
            fP1_DP = sum(Ref_DP_final)/(sum(Sum_DP_final)))

# final depths after removing dubious Het calls
final_allele_depth <- SNPs_merge_finalGT %>% 
  filter(!Cut, GT == "0/1", GQ >= 50, !(GT == "0/1" & Alt_DP_final <= 3)) %>% 
  ungroup() %>% group_by(ID) %>% 
  summarize(call = "Final", 
            tot_P1_DP = sum(Ref_DP_final), 
            tot_P2_DP = sum(Alt_DP_final),
            fP1_DP = sum(Ref_DP_final)/(sum(Sum_DP_final)))

allele_depth <- rbind(P1call_allele_depth, P2call_allele_depth, merge_allele_depth, final_allele_depth)
allele_depth <- allele_depth %>% mutate(total_DP = (tot_P1_DP + tot_P2_DP))
allele_depth$call <- factor(allele_depth$call, levels = c("P1", "P2", "Merge", "Final"))

fP1_P1cVsP2c_perm <- perm_test(allele_depth, cat_var = "call", cat_names = c("P1", "P2"), 
                               response_var = "fP1_DP", test_stat = mean, n_perms = 10000)
fP1_P1cVsP2c_perm <- c(fP1_P1cVsP2c_perm, cat_1 = "P1", cat_2 = "P2")

fP1_P1cVsMerge_perm <- perm_test(allele_depth, cat_var = "call", cat_names = c("P1", "Merge"), 
                                 response_var = "fP1_DP", test_stat = mean, n_perms = 10000)
fP1_P1cVsMerge_perm <- c(fP1_P1cVsMerge_perm, cat_1 = "P1", cat_2 = "Merge")

fP1_P2cVsMerge_perm <- perm_test(allele_depth, cat_var = "call", cat_names = c("P2", "Merge"), 
                                 response_var = "fP1_DP", test_stat = mean, n_perms = 10000)
fP1_P2cVsMerge_perm <- c(fP1_P2cVsMerge_perm, cat_1 = "P2", cat_2 = "Merge")

fP1_MergeVsfinal_perm <- perm_test(allele_depth, cat_var = "call", cat_names = c("Merge", "Final"), 
                                   response_var = "fP1_DP", test_stat = mean, n_perms = 10000)
fP1_MergeVsfinal_perm <- c(fP1_MergeVsfinal_perm, cat_1 = "Merge", cat_2 = "Final")

fP1_perm_merge <- as.data.frame(bind_rows(fP1_P1cVsP2c_perm, fP1_P1cVsMerge_perm, 
                                      fP1_P2cVsMerge_perm, fP1_MergeVsfinal_perm))
fP1_perm_merge$cat_1 <- factor(fP1_perm_merge$cat_1, levels = c("P1", "P2", "Merge", "Final"))
fP1_perm_merge$cat_2 <- factor(fP1_perm_merge$cat_2, levels = c("P1", "P2", "Merge", "Final"))

fP1_perm_merge <- cbind(fP1_perm_merge, sig = ifelse(fP1_perm_merge$rejectNull, "*", " "))
# allele_depth %>% filter(call == "P1") %>% pull("fP1_DP") %>% t.test(., mu = 0.5)
# allele_depth %>% filter(call == "P2") %>% pull("fP1_DP") %>% t.test(., mu = 0.5)
# allele_depth %>% filter(call == "Merge") %>% pull("fP1_DP") %>% t.test(., mu = 0.5)
# allele_depth %>% filter(call == "Final") %>% pull("fP1_DP") %>% t.test(., mu = 0.5)

# Plot frequency of P1 reference allele among callsets colored by total depth
###############################################################################
DP_mid <- median(allele_depth$total_DP)
y_range <- 0.04
allele_mapping_bias_plot <- allele_depth %>% 
  ggplot() + 
  geom_hline(aes(yintercept = 0.5), size = 0.2) +

  geom_jitter(aes(x = call, y = fP1_DP, color = total_DP), 
              height = 0, width = 0.4, alpha = 0.7, size = 2) +
  geom_boxplot(aes(x = call, y = fP1_DP), outlier.colour = NA, width = 0.5, fill = NA) +
  # geom_bracket(data = fP1_perm_merge %>% filter(rejectNull), aes(xmin = cat_1, xmax = cat_2, label = sig), 
  #              y.position = 0.5 - 0.01 + y_range, step.increase = 0.075, vjust = 0.7, label.size = 5, fontface = "bold") +
  ylim(0.5 - y_range, 0.5 + y_range) +
  scale_color_gradient2(low = "orangered4", mid = "grey70", high = "blue1", midpoint = DP_mid + 1.5E6) +
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
  summarize(fRef_P1 = Ref_DP_P1call/(Ref_DP_P1call + Alt_DP_P1call),
            fRef_P2 = Ref_DP_P2call/(Ref_DP_P2call + Alt_DP_P2call),
            fRef_final = Ref_DP_final/(Ref_DP_final + Alt_DP_final),
            het_cut = Cut) 

allele_freq %>%
  ggplot() + 
  geom_histogram(aes(x = fRef_P1), fill = "brown2", color = "brown2",
                 binwidth = 0.025, alpha = 0.5) +
  geom_histogram(aes(x = fRef_P2), fill = "dodgerblue2", color = "dodgerblue2",
                 binwidth = 0.025, alpha = 0.5) +
  geom_histogram(aes(x = fRef_final), fill = "white", color = "grey20",
                 binwidth = 0.025, alpha = 0.3)

  
allele_depth %>% 
  ggplot() + geom_abline(color = "blue4") +
  geom_point(aes(x = tot_P1_DP/1E6, y = tot_P2_DP/1E6, shape = call, group = ID), size = 2) + 
  geom_line(aes(x = tot_P1_DP/1E6, y = tot_P2_DP/1E6, group = ID), color = "grey30", alpha = 0.5) +
  # geom_label_repel(data = subset(allele_depth, tot_P1_DP/tot_P2_DP < 0.8 & call == "P2"), 
  #                  aes(x = tot_P1_DP, y = tot_P2_DP, label = ID), label.size = 0, force = 1.3) + 
  xlab("Reads supporting P1 allele") + ylab("Reads supporting P2 allele") +
  scale_shape_manual(values = c(15, 16, 17, 18), name = "Reference")

allele_depth %>% 
  ggplot() + geom_histogram(aes(x = fP1_DP)) + facet_wrap(~call, ncol = 1)

SNPs_merge_finalGT %>% 
  filter(GT == "0/1", !Cut) %>%
  ggplot() + geom_histogram(aes(x = Ref_DP_final/(Ref_DP_final + Alt_DP_final)),
                            breaks = seq(0.025, 0.975, 0.05)) + 
  facet_wrap(~Line)

# Compare distribution of P1 allele frequencies for het sites and dubious het sites
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
  select(CHROM:POSi, ID, GT_P1call, GQ_P1call, GT_P2call, GQ_P2call, GT, GQ) %>% View()
  dplyr::count(GT, GT_P1call, GT_RMcall)

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
