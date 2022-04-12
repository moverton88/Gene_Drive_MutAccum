

###############################################################################
# Mean LOH rates and permutation tests for each chromosome ------
# Calculated as n_LOH_events/n_clones/gens
POSi_data_in <- all_LOHbounds_merge %>% 
  filter(length > 1) # %>% 
  # filter(!isTerm) %>%
  # filter(!(start_POSi %in% prob_sites | end_POSi %in% prob_sites)) %>%
  # select(-is_error)

POSi_het <- all_GT_bounds_merge %>% 
  filter(GT == "0/1") # %>% 
  # filter(!(start_POSi %in% prob_sites | end_POSi %in% prob_sites)) # %>%
  # select(-is_error)

POSi_data_in$est_start_POS <- ConvertPosIndicies(POSi_data_in, pos_col = "est_start", index_out = "POS")
POSi_data_in$est_end_POS <- ConvertPosIndicies(POSi_data_in, pos_col = "est_end", index_out = "POS")

POSi_data_in$dist_cent_start <- ChromosomeCoordinates(POSi_data_in, POSi_col = "est_start") %>% 
  pull(dist_cent)
POSi_data_in$dist_cent_end <- ChromosomeCoordinates(POSi_data_in, POSi_col = "est_end") %>% 
  pull(dist_cent)
POSi_data_in <- POSi_data_in %>% mutate(est_mid_POS = (est_start_POS + est_end_POS)/2)

POSi_data_in$ID_CHROM <- paste0(POSi_data_in$ID, "_", POSi_data_in$CHROM)
whl_chrom <- POSi_data_in %>% 
  filter(dist_cent_start < 0, dist_cent_end > 0, isTerm) %>% 
  select(ID, CHROM, GT)

wh_chrom_LOH <- all_LOHbounds_merge_EC %>% 
  count(ID, CHROM, GT) %>% 
  pivot_wider(names_from = GT, values_from = n, values_fill = 0)

colnames(wh_chrom_LOH)[3:5] <- c("het", "hom_ref", "hom_alt")
wh_chrom_ID_CHROM <- wh_chrom_LOH %>% 
  filter(het == 0 & (hom_ref == 1 | hom_alt == 1)) %>% 
  mutate(ID_CHROM = paste0(ID, "_", CHROM)) %>% pull(ID_CHROM)
i_wh_chrom <- POSi_data_in$ID_CHROM %in% wh_chrom_ID_CHROM

i_Term <- POSi_data_in$isTerm
head_term <- POSi_data_in$dist_cent_start < 0 & POSi_data_in$dist_cent_end < 0
tail_term <- POSi_data_in$dist_cent_start > 0 & POSi_data_in$dist_cent_end > 0
cross_h_term <- POSi_data_in$dist_cent_start < 0 & POSi_data_in$dist_cent_end > 0
cross_ID_CHROM <- POSi_data_in[i_Term & cross_h_term & !i_wh_chrom, ] %>% pull(ID_CHROM)

cross_h_term_ID <- POSi_data_in %>% 
  filter(isTerm) %>%
  # mutate(ID_CHROM = paste0(ID, "_", CHROM)) %>% 
  filter(ID_CHROM %in% cross_ID_CHROM) %>% distinct(ID_CHROM, .keep_all = T) %>% 
  summarize(ID_CHROM = ID_CHROM, head_term = ifelse(est_start_POS == 1, T, F))

# Correct POS for terminal LOH that do not cross the centromere
POSi_data_in$est_mid_POS[i_Term] <- ifelse(head_term[i_Term], 
                                           POSi_data_in$est_end_POS[i_Term], 
                                           POSi_data_in$est_start_POS[i_Term])

# Correct POS for terminal LOH that do cross the centromere
i_crct <- i_Term & cross_h_term & !i_wh_chrom
POSi_data_in$est_mid_POS[i_crct][cross_h_term_ID$head_term] <- 
  POSi_data_in$est_end_POS[i_crct][cross_h_term_ID$head_term]
POSi_data_in$est_mid_POS[i_crct][!cross_h_term_ID$head_term] <- 
  POSi_data_in$est_start_POS[i_crct][!cross_h_term_ID$head_term]
POSi_data_in$est_mid <- ConvertPosIndicies(POSi_data_in, pos_col = "est_mid_POS", index_out = "POSi")

POSi_data_wHet <- merge(POSi_data_in, POSi_het, all = T) %>% arrange(ID, est_start)
# 
# POSi_data_wHet <- all_LOHbounds_2_merge_EC %>% 
#   filter(!is_error) %>% select(-is_error)
# 
# POSi_data_wHet$est_start_POS <- ConvertPosIndicies(POSi_data_wHet, 
#                                                    pos_col = "est_start", index_out = "POS")
# POSi_data_wHet$est_end_POS <- ConvertPosIndicies(POSi_data_wHet, 
#                                                  pos_col = "est_end", index_out = "POS")
# POSi_data_wHet$est_mid <- ConvertPosIndicies(POSi_data_wHet, 
#                                              pos_col = "est_mid_POS", index_out = "POSi")

# POSi_data_wHet$est_start_POS[i_crct] <- 1

LOH_bounds_plot <- POSi_data_wHet %>%
# LOH_bounds_plot <- POSi_data_in %>%
  # all_LOHbounds_2_merge_EC %>%
  # filter(!is_error) %>% 
  # filter(!(start_POSi %in% prob_sites | end_POSi %in% prob_sites)) %>%
  # filter(GT != "0/1") %>%
  # filter(est_start > 5800000) %>%
  filter(CHROM == "04") %>%
  # filter(Tx == "F") %>%
  ggplot() + 
  geom_segment(aes(x = est_start, xend = est_end, 
                   y = ID, yend = ID, color = GT), size = 0.75) +
  # geom_point(aes(x = est_mid, y = ID, color = GT), size = 0.5) +
  scale_color_manual(values = c(allelePal[1], "grey80", allelePal[3])) +
  # scale_color_manual(values = allelePal[c(1,3)]) +
  scale_x_continuous(expand = c(0,0)) +
  facet_wrap(~CHROM, scales = "free_x")

LOH_bounds_plot

ggsave(file.path(outIntDir, "LOH_bounds_plot.png"), 
       plot = LOH_bounds_plot,
       device = "png",
       width = 50, height = 30, 
       units = "in",
       dpi = 600, limitsize = FALSE)

# Use the LOHrateSW() function to determine clone/chromosomes too few markers
all_LOHrates_Chrm <- POSi_data_wHet %>% 
  LOHrateSW(by_GT = F, win = "CHROM", 
            chrom_df = chrom_bound_BY,
            rm_zero_valids = F)

low_cover_ID_CHROM <- all_LOHrates_Chrm %>% 
  filter(prop_valid < 0.4) %>% select(ID, CHROM) %>% 
  mutate(ID_CHROM = paste0(ID, "_", CHROM))

POSi_data_in$ID_CHROM <- paste0(POSi_data_in$ID, "_", POSi_data_in$CHROM)

POSi_data_in$lc_CHROM <- POSi_data_in$ID_CHROM %in% low_cover_ID_CHROM$ID_CHROM

# Calculate LOH event rates
LOHcounts_ID_CHROM <- POSi_data_in %>% 
  filter(!lc_CHROM) %>%
  group_by(ID, CHROM, .drop = F) %>% 
    count(name = "n_LOH") %>% as.data.frame()

LOHcounts_ID_CHROM <- CategoriesFromID(LOHcounts_ID_CHROM)
LOHcounts_ID_CHROM$Tx_ID <- Recode_Tx_ID(LOHcounts_ID_CHROM$Tx)

LOHcounts_CHROM_mean <- LOHcounts_ID_CHROM %>% 
  group_by(CHROM) %>% 
  summarize(total_LOH = sum(n_LOH), mean_LOH = mean(n_LOH), se_LOH = se(n_LOH), 
            mean_rate = mean(n_LOH)/n_gens, se_rate = se(n_LOH)/n_gens, 
            CI_rate = 1.96*se(n_LOH)/n_gens)

LOHcounts_CHROM_mean$CHROM_len <- chrom_lengths_BY
LOHcounts_CHROM_mean <- LOHcounts_CHROM_mean %>% mutate(bp_rate = total_LOH/CHROM_len)

max(LOHcounts_CHROM_mean$bp_rate)/min(LOHcounts_CHROM_mean$bp_rate)

LOHcounts_ID_CHROM_mean <- LOHcounts_ID_CHROM %>% 
  group_by(Tx_ID, CHROM) %>% 
  summarize(total_LOH = sum(n_LOH), mean_LOH = mean(n_LOH), se_LOH = se(n_LOH), 
            mean_rate = mean(n_LOH)/n_gens, se_rate = se(n_LOH)/n_gens, 
            CI_rate = 1.96*se(n_LOH)/n_gens)

LOHcounts_ID_CHROM_mean$rom_CHROM <- roman_chr[as.numeric(LOHcounts_ID_CHROM_mean$CHROM)]

Tx_combo <- data.frame(Tx_1 = c("WT", "WT", "Cas9"), 
                       Tx_2 = c("Cas9", "Drive", "Drive"))

all_CHROM_LOHrate_perm <- data.frame(NULL)
for(tx in 1:nrow(Tx_combo)){
  # tx = 1
  all_Chr_perm <- data.frame(Combo = with(Tx_combo[tx, ], paste(Tx_1, Tx_2, sep = "-")), 
                             CHROM = chrom_bound_BY$CHROM, 
                             obsvStat = 0, critVal = 0, pVal = 0, rejectNull = 0)
  
  for(c in seq_along(chrom_bound_BY$CHROM)) {
    # c = 14
    chr = chrom_bound_BY$CHROM[c]
    Chr_LOHcounts_ID <- LOHcounts_ID_CHROM %>% filter(CHROM == chr)
    Chr_perm <- perm_test(Chr_LOHcounts_ID, 
                          cat_var = "Tx_name", 
                          cat_names = with(Tx_combo[tx, ], c(Tx_1, Tx_2)), 
                          response_var = "n_LOH", 
                          alpha = 0.05,
                          n_perms = 10000, alt_hyp = "two-tailed",
                          rtrn = "all", include_matrix = F)
    all_Chr_perm[c, 3:6] <- unlist(Chr_perm[1:4])
  }
  all_CHROM_LOHrate_perm <- rbind(all_CHROM_LOHrate_perm, all_Chr_perm)
}

# all_CHROM_LOHrate_perm %>% arrange(desc(obsvStat)) %>% head()
# all_CHROM_LOHrate_perm %>% arrange(pVal) %>% head()

all_CHROM_LOHrate_perm_slim <- all_CHROM_LOHrate_perm %>% filter(Combo != "Cas9-Drive")
all_CHROM_LOHrate_perm_slim <- BHcorrection(all_CHROM_LOHrate_perm_slim)

# all_CHROM_LOHrate_perm <- BHcorrection(all_CHROM_LOHrate_perm)


CHROM_perm_sig <- all_CHROM_LOHrate_perm %>% filter(pVal < FDR)
all_CHROM_LOHrate_perm$sig_lab_p <- ifelse(all_CHROM_LOHrate_perm$rejectNull == 1, "*", "")
all_CHROM_LOHrate_perm$sig_lab_BH <- ifelse(all_CHROM_LOHrate_perm$rejectNull_BH == 1, "**", "")
all_CHROM_LOHrate_perm$sig_lab_both <- paste0(all_CHROM_LOHrate_perm$sig_lab_p, 
                                                ifelse(all_CHROM_LOHrate_perm$rejectNull_BH == 1, ", ", ""),
                                              all_CHROM_LOHrate_perm$sig_lab_BH)

combo_cols <- colsplit(all_CHROM_LOHrate_perm$Combo, "-", c("Tx1", "Tx2"))
all_CHROM_LOHrate_perm <- cbind(all_CHROM_LOHrate_perm, combo_cols)

# all_CHROMrate_perm_order %>% ggplot(aes(x = pVal)) + geom_histogram(binwidth = 0.05)

all_CHROM_LOHrate_perm$Tx1 <- factor(all_CHROM_LOHrate_perm$Tx1, levels = c("WT", "Cas9", "Drive"))
all_CHROM_LOHrate_perm$Tx2 <- factor(all_CHROM_LOHrate_perm$Tx2, levels = c("WT", "Cas9", "Drive"))
all_CHROM_LOHrate_perm$Combo <- factor(all_CHROM_LOHrate_perm$Combo)
sig_CHROMrate <- all_CHROM_LOHrate_perm %>% filter(rejectNull_BH == 1) 

LOH_xCHROM_aov <- aov(n_LOH ~ Tx_name + CHROM, data = LOHcounts_ID_CHROM)
summary(LOH_xCHROM_aov)
LOH_xCHROM_Tuk <- TukeyHSD(LOH_xCHROM_aov, "Tx_name")
LOH_xCHROM_Tuk

y_loc <- max(LOHcounts_ID_CHROM_mean$mean_rate + LOHcounts_ID_CHROM_mean$CI_rate)
LOHrate_CHROM_line <- LOHcounts_ID_CHROM_mean %>% 
  # filter(Tx_name != "Cas9") %>%
  ggplot() + 
  # geom_errorbar(aes(x = Tx_name, ymin = mean_rate, ymax = mean_rate, color = Tx_name), 
  #               width = 0.6, position = position_dodge(width = 0.6), size = 0.75) +
  geom_errorbar(aes(x = Tx_ID, ymin = ifelse(mean_rate - CI_rate < 0, 0, mean_rate - CI_rate), 
                    ymax = mean_rate + CI_rate, color = Tx_ID),
                width = 0.2, position = position_dodge(width = 0.6), size = 0.75) +
  geom_point(aes(x = Tx_ID, y = mean_rate, color = Tx_ID), 
                position = position_dodge(width = 0.6), size = 4) +
  # geom_bracket(data = sig_CHROMrate,
               # y.position = y_loc*1.05 + (as.numeric(sig_CHROMrate$Combo)-1) * y_loc * 0.03, 
                # step.increase = 0,
  #              vjust = 0.4, label.size = 5,
  #              aes(xmin = Tx1, xmax = Tx2, label = sig_lab_BH)) +
  # geom_text(aes(x = Tx_name, y = -0.00005, label = round(mean_LOH, 2), color = Tx_name),
  #           position = position_dodge(width = 0.6), show.legend = FALSE) +
  scale_color_manual(values = txPal[c(3, 2, 1)], name = "Strain") +
  scale_y_continuous(labels = c(0, format(seq(2.5e-4, 1e-3, 2.5e-4), scientific = T)),
                     breaks = c(0, seq(2.5e-4, 1e-3, 2.5e-4))) +
  # ylim(c(0, NA)) +
  ylab("LOH event rate (/clone/generation)") + 
  xlab("Chromosome") +
  facet_wrap(~rom_CHROM, ncol = 8, strip.position = "bottom") +
  theme(axis.text.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.background = element_rect(color = "grey80"),
        text = element_text(size = 18), 
        legend.position = c(0.03, 0.93),
        legend.background = element_rect(fill = "white", color = "white"))

LOHrate_CHROM_line

LOHcounts_ID_CHROM_mean$n_CHROM <- as.numeric(LOHcounts_ID_CHROM_mean$CHROM)
LOHcounts_ID_CHROM_mean$Tx_n <- as.numeric(LOHcounts_ID_CHROM_mean$Tx_ID)

point_diff <- 0.16

LOHrate_CHROM_line <- LOHcounts_ID_CHROM_mean %>% 
  ggplot() + 
  annotate(geom = "rect", xmin = LOHcounts_ID_CHROM_mean$n_CHROM[c(F, T)] - 0.5, 
           xmax = LOHcounts_ID_CHROM_mean$n_CHROM[c(F, T)] + 0.5, ymin = 0, ymax = Inf, 
           fill = "grey90", alpha = 0.1) +
  # geom_errorbar(aes(x = Tx_name, ymin = mean_rate, ymax = mean_rate, color = Tx_name), 
  #               width = 0.6, position = position_dodge(width = 0.6), size = 0.75) +
  geom_line(aes(x = n_CHROM + (Tx_n - 2) * point_diff, y = mean_rate, color = Tx_ID, group = Tx_ID), 
             alpha = 0.7, size = 1) +
  geom_errorbar(aes(x = n_CHROM + (Tx_n - 2) * point_diff, 
                    ymin = ifelse(mean_rate - CI_rate < 0, 0, mean_rate - CI_rate), 
                    ymax = mean_rate + CI_rate, color = Tx_ID),
                width = 0, alpha = 0.7, size = 0.5) +
  geom_point(aes(x = n_CHROM + (Tx_n - 2) * point_diff, y = mean_rate, color = Tx_ID), 
             size = 4) +
  # geom_bracket(data = sig_CHROMrate,
  #              y.position = y_loc*1.05 + (as.numeric(sig_CHROMrate$Combo)-1) * y_loc * 0.03, step.increase = 0,
  #              vjust = 0.4, label.size = 5,
  #              aes(xmin = Tx1, xmax = Tx2, label = sig_lab_BH)) +
  # geom_text(aes(x = Tx_name, y = -0.00005, label = round(mean_LOH, 2), color = Tx_name),
  #           position = position_dodge(width = 0.6), show.legend = FALSE) +
  scale_color_manual(values = txPal[c(3, 2, 1)], name = "Strain") +
  scale_y_continuous(labels = c(0, format(seq(2.5e-4, 1e-3, 2.5e-4), scientific = T)),
                     breaks = c(0, seq(2.5e-4, 1e-3, 2.5e-4))) +
  scale_x_continuous(labels = roman_chr, expand = c(0, 0),
                     breaks = 1:16, limits = c(0.5, 16.5)) +
  # ylim(c(0, NA)) +
  ylab("LOH event rate (/clone/generation)") + 
  xlab("Chromosome") +
  # facet_wrap(~rom_CHROM, ncol = 8, strip.position = "bottom") +
  theme(axis.text.x = element_text(vjust = 5),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        # panel.background = element_rect(color = "grey80"),
        text = element_text(size = 24),
        legend.position = c(0.075, 0.87),
        legend.background = element_rect(fill = "white", color = "white"))

LOHrate_CHROM_line

ggsave(file.path(outIntDir, "LOHrate_CHROM_line_2022_03.png"), 
       plot = LOHrate_CHROM_line,
       device = "png",
       width = 16, height = 9, 
       units = "in",
       dpi = 600)

###############################################################################
# Empirical cumulative distribution of LOH rates along the length of each 
# chromosome with permuted Anderson-Darling test
POSi_test_in <- POSi_data_in %>%
  filter(!isTerm)

all_CHROM_LOHprop_KS <- data.frame(NULL)
for(tx in 1:nrow(Tx_combo)){
  # tx = 1
  Tx_1 <- Tx_combo$Tx_1[tx]
  Tx_2 <- Tx_combo$Tx_2[tx]
  Combo_Chr_KS <- data.frame(Combo = with(Tx_combo[tx, ], paste(Tx_1, Tx_2, sep = "-")),
                             Tx_1 = Tx_1, Tx_2 = Tx_2,
                             CHROM = chrom_bound_BY$CHROM, 
                             KSstat = 0, pVal = 0)
  KS_LOHcounts_1 <- POSi_test_in %>% filter(Tx_name %in% Tx_1) %>% 
    select(Tx_name, ID, CHROM, est_mid_POS)
  KS_LOHcounts_2 <- POSi_test_in %>% filter(Tx_name %in% Tx_2) %>% 
    select(Tx_name, ID, CHROM, est_mid_POS)
  for(ch in seq_along(chrom_bound_BY$CHROM)) {
    # ch = 1
    chr = chrom_bound_BY$CHROM[ch]
    Chr_LOHcounts_1 <- KS_LOHcounts_1 %>% filter(CHROM == chr)
    Chr_LOHcounts_2 <- KS_LOHcounts_2 %>% filter(CHROM == chr)
    Chr_KS <- ks.test(Chr_LOHcounts_1$est_mid_POS, Chr_LOHcounts_2$est_mid_POS, alternative = "two.sided")
    Combo_Chr_KS[ch, c(5, 6)] <- c(Chr_KS$statistic, Chr_KS$p.value)
  }
  all_CHROM_LOHprop_KS <- rbind(all_CHROM_LOHprop_KS, Combo_Chr_KS)
}

# all_CHROM_LOHprop_KS <- BHcorrection(all_CHROM_LOHprop_KS)
all_CHROM_LOHprop_KS_slim <- all_CHROM_LOHprop_KS %>% filter(Combo != "Cas9-Drive")
all_CHROM_LOHprop_KS_slim <- BHcorrection(all_CHROM_LOHprop_KS_slim)

all_CHROM_LOHprop_AD <- data.frame(NULL)
for(tx in 1:nrow(Tx_combo)){
  # tx = 1
  Tx_1 <- Tx_combo$Tx_1[tx]
  Tx_2 <- Tx_combo$Tx_2[tx]
  Combo_Chr_AD <- data.frame(Combo = with(Tx_combo[tx, ], paste(Tx_1, Tx_2, sep = "-")),
                             Tx_1 = Tx_1, Tx_2 = Tx_2,
                             CHROM = chrom_bound_BY$CHROM, 
                             ADstat = 0, pVal = 0)
  AD_LOHcounts <- POSi_test_in %>% filter(Tx_name %in% c(Tx_1, Tx_2)) %>% 
    # filter(!isTerm) %>%
    select(Tx_name, ID, CHROM, est_mid_POS)
  AD_LOHcounts$Tx_name <- droplevels(AD_LOHcounts$Tx_name)
  for(ch in seq_along(chrom_bound_BY$CHROM)) {
    # ch = 1
    chr = chrom_bound_BY$CHROM[ch]
    Chr_LOHcounts_ID <- AD_LOHcounts %>% filter(CHROM == chr)
    Chr_AD <- ad.test(est_mid_POS ~ Tx_name, data = Chr_LOHcounts_ID, method = "exact")
    Combo_Chr_AD[ch, c(5, 6)] <- Chr_AD$ad[2, c(1, 3)]
  }
  all_CHROM_LOHprop_AD <- rbind(all_CHROM_LOHprop_AD, Combo_Chr_AD)
}

all_CHROM_LOHprop_AD_slim <- all_CHROM_LOHprop_AD %>% filter(Combo != "Cas9-Drive")
all_CHROM_LOHprop_AD_slim <- BHcorrection(all_CHROM_LOHprop_AD_slim)
ad_label <- all_CHROM_LOHprop_AD_slim %>% filter(rejectNull_BH == 1)
ad_label$sig <- "*"
ad_label$Tx_2 <- factor(ad_label$Tx_2, levels = Tx_name_levels)
ad_label$rom_CHROM <- factor(chrom_IDs$rom_CHROM[as.numeric(ad_label$CHROM)],
                             levels = levels(chrom_IDs$rom_CHROM))
nLOH_CHROM <- POSi_data_in %>% count(Tx_name, CHROM)

marker_set <- data.frame(POSi = clean_markers)
marker_set_POS <- ConvertPosIndicies(marker_set, pos_col = "POSi", 
                                     index_out = "POS", add_chroms = T)
marker_set <- cbind(marker_set, marker_set_POS)
POSi_data_in$rom_CHROM <- factor(chrom_IDs$rom_CHROM[as.numeric(POSi_data_in$CHROM)],
                                 levels = levels(chrom_IDs$rom_CHROM))

# ad_label gives chromosomes and treatment pairs that are significant
# get position information to make brackets for indicating significance on plot

ecd_sig_region <- POSi_data_in %>%
  filter(!isTerm, rom_CHROM %in% ad_label$rom_CHROM, Tx_name %in% c("WT", "Drive")) %>% 
  arrange(Tx_name, est_mid_POS) %>% group_by(Tx_name) %>% 
  summarize(rom_CHROM = rom_CHROM, POS = est_mid_POS, f_bp = row_number()/n()) 

# Chr V is significant b/t WT and Drive. Cannot program bracket values
# for more than one significant combination, so this is coded for only one.
# find position that is maximally different between strains
ecd_sig_region <- ecd_sig_region %>% ungroup() %>% arrange(POS) %>% 
  mutate(lag_diff = f_bp - lag(f_bp))

ecd_sig_hi_POS <- ecd_sig_region %>% arrange(desc(lag_diff)) %>% head(1)
max_POS <- ecd_sig_hi_POS$POS + 5000

ecd_sig_region_lo <- ecd_sig_region %>% filter(Tx_name != ecd_sig_hi_POS$Tx_name)
# Find corresponding POS for lower value treatment
ecd_sig_lo_POS <- ecd_sig_region_lo[which.min(abs(max_POS - ecd_sig_region_lo$POS)), ]
ecd_bracket <- data.frame(rom_CHROM = ecd_sig_hi_POS$rom_CHROM, 
                          POS = b_POS, y_lo = ecd_sig_lo_POS$f_bp, y_hi = ecd_sig_hi_POS$f_bp)
ecd_bracket <- rbind(ecd_sig_hi_POS, ecd_sig_lo_POS)
ecd_bracket$POS <- mean(ecd_bracket$POS)

LOH_CHROM_ecd <- POSi_data_in %>%
  filter(!isTerm) %>%
  ggplot() +
  geom_segment(data = chrom_lengths_BY_df, 
               aes(x = 0, y = 0, xend = chrom_length/1000, yend = 1), size = 0.1, alpha = 0.9) +
  # geom_polygon(data = ecd_sig_region, aes(x = POS/1000, y = n, fill = Tx_name)) +
  stat_ecdf(aes(x = est_mid_POS/1000, color = Tx_name)) + 
  # geom_segment(data = ecd_bracket,
  #              aes(x = POS/1000, xend = POS/1000, y = y_lo, yend = y_hi)) +
  # geom_label(data = ad_label, aes(x = 10 , y = 0.85, label = sig, color = Tx_2), 
  #            size = 9, hjust = "left", label.size = 0, show.legend = F) +
  # geom_point(data = marker_set, aes(x = POS/1000, y = -0.05), shape = "|", size = 0.5) +
  # stat_ecdf(aes(x = POSi), color = "red", size = 0.3) + 
  scale_color_manual(values = txPal, name = "Strain") + 
  xlab("Chromosome Position (kbp)") + ylab("Cumulative Fraction of LOH events") +
  scale_x_continuous(expand = c(0.01, 0.01)) +
  facet_wrap(~rom_CHROM, ncol = 4, scales = "free_x") +
  theme(
    panel.grid.minor.y = element_blank(), 
    panel.grid.minor.x = element_blank(),
    panel.background = element_rect(color="grey80"),
    text = element_text(size = 18), legend.position = "bottom")

LOH_CHROM_ecd

ggsave(file.path(outIntDir, "iLOH_BPrate_xCHROM_ecd_2022_02.png"), 
       plot = LOH_CHROM_ecd,
       device = "png",
       width = 16, height = 9, 
       units = "in",
       dpi = 600)

chrom_in <- c("02", "03", "05", "09")

ecd_in <- POSi_data_in %>% filter(CHROM %in% chrom_in)
ecd_in$CHROM <- factor(ecd_in$CHROM, levels = chrom_in)
marker_in <- marker_set %>% filter(CHROM %in% chrom_in)
marker_in$CHROM <- factor(marker_in$CHROM, levels = chrom_in)
chrom_lengths_in <- chrom_lengths_BY_df %>% filter(CHROM %in% chrom_in)
chrom_lengths_in$CHROM <- factor(chrom_lengths_in$CHROM, levels = chrom_in)

LOH_CHROM_ecd_sub <- ecd_in %>%
  ggplot() +
   # geom_label(aes(x = chrom_bound_BY$End[1], y = 0.85, label = ad_label), hjust = "left") +
  # facet_wrap(~CHROM, ncol = 2, scales = "free_x") +
  stat_ecdf(aes(x = est_mid_POS/1000, color = Tx_name)) + 
  geom_point(data = marker_in, aes(x = POS/1000, y = -0.05), shape = "|", size = 1) +
  geom_segment(data = chrom_lengths_in, 
               aes(x = 0, y = 0, xend = chrom_length/1000, yend = 1), size = 0.1, alpha = 0.9) +
  # stat_ecdf(aes(x = POSi), color = "red", size = 0.3) + 
  scale_color_manual(values = txPal, name = "Strain") + 
  xlab("Chromosome Position (kb)") + ylab("Cumulative Fraction of Breakpoints") +
  scale_x_continuous(expand = c(0.01, 0.01)) +
  facet_wrap(~CHROM, ncol = 2, scales = "free_x") +
  theme(
        panel.grid.minor.y = element_blank(), 
        panel.grid.minor.x = element_blank(),
        panel.background = element_rect(color="grey80"),
        text = element_text(size = 16), legend.position = "bottom")

LOH_CHROM_ecd_sub

ggsave(file.path(outIntDir, "iLOH_BPrate_ecd_sub_2022_02.png"), 
       plot = LOH_CHROM_ecd_sub,
       device = "png",
       width = 11, height = 8.5, 
       units = "in",
       dpi = 600)


# Fractional distribution of LOH events by chromosome #########################
POSi_test_in <- POSi_data_in %>%
  filter(!isTerm)

POSi_test_in


###############################################################################
# LOH event rate for 1/3 chromosomes

all_LOHrates_Chrm3 <- POSi_data_wHet %>% LOHrateSW(by_GT = F, win = "CHROM/3", chrom_df = chrom_bound_BY)

low_cover_ID_CHROM3 <- all_LOHrates_Chrm3 %>% 
  filter(prop_valid < 0.3) %>% select(ID, CHROM3) %>% 
  mutate(ID_CHROM3 = paste0(ID, "_", CHROM3))

chrom3s <- floor(chrom_lengths_BY/3)

chrom3_bound_BY <- data.frame(CHROM = chrom_bound_BY$CHROM, 
                              Start = chrom_bound_BY$Start,
                              mid1 = chrom_bound_BY$Start + chrom3s,
                              mid2 = chrom_bound_BY$Start + chrom3s*2,
                              End = chrom_bound_BY$End)

POSi_data_in$CHROM3 <- "00"
for(i in 1:16) {
  # i = 1
  chr <- chrom3_bound_BY[i, "CHROM"]
  chr3 <- paste(chr, 1:3, sep = "_")
  i_pos1 <- POSi_data_in$est_mid >= chrom3_bound_BY[i, "Start"] & 
    POSi_data_in$est_mid < chrom3_bound_BY[i, "mid1"]
  i_pos2 <- POSi_data_in$est_mid >= chrom3_bound_BY[i, "mid1"] & 
    POSi_data_in$est_mid < chrom3_bound_BY[i, "mid2"]
  i_pos3 <- POSi_data_in$est_mid >= chrom3_bound_BY[i, "mid2"] & 
    POSi_data_in$est_mid <= chrom3_bound_BY[i, "End"]
  POSi_data_in$CHROM3[i_pos1] <- chr3[1]
  POSi_data_in$CHROM3[i_pos2] <- chr3[2]
  POSi_data_in$CHROM3[i_pos3] <- chr3[3]
}

POSi_data_in$CHROM3 <- factor(POSi_data_in$CHROM3)
POSi_data_in$ID_CHROM3 <- paste0(POSi_data_in$ID, "_", POSi_data_in$CHROM3)

POSi_data_in$lc_CHROM3 <- POSi_data_in$ID_CHROM3 %in% low_cover_ID_CHROM3$ID_CHROM3

LOHcounts_ID_CHROM3 <- POSi_data_in %>% 
  filter(!lc_CHROM3) %>%
  group_by(ID, CHROM3, .drop = F) %>% count(name = "n_LOH") %>% as.data.frame()
LOHcounts_ID_CHROM3 <- CategoriesFromID(LOHcounts_ID_CHROM3)
LOHcounts_ID_CHROM3_mean <- LOHcounts_ID_CHROM3 %>% 
  group_by(Tx_name, CHROM3) %>% 
  summarize(total_LOH = sum(n_LOH), mean_LOH = mean(n_LOH), se_LOH = se(n_LOH), 
            mean_rate = mean(n_LOH)/n_gens, se_rate = se(n_LOH)/n_gens, CI_rate = 1.96*se(n_LOH)/n_gens)
LOHcounts_ID_CHROM3_mean$CHROM <- factor(substr(LOHcounts_ID_CHROM3_mean$CHROM3, 1, 2))
LOHcounts_ID_CHROM3_mean$CHROM3_pos <- factor(substr(LOHcounts_ID_CHROM3_mean$CHROM3, 4, 4))

ch3_names <- sort(unique(LOHcounts_ID_CHROM3$CHROM3))
all_CHROM3_LOHrate_perm <- data.frame(NULL)
for(tx in 1:nrow(Tx_combo)){
  # tx = 1
  all_chr3_perm <- data.frame(Combo = with(Tx_combo[tx, ], paste(Tx_1, Tx_2, sep = "-")), 
                              CHROM3 = ch3_names, 
                             obsvStat = 0, critVal = 0, pVal = 0, rejectNull = 0)
  for(c_i in seq_along(ch3_names)) {
    # chr = 10
    chr3 = ch3_names[c_i]
    chr3_LOHcounts_ID <- LOHcounts_ID_CHROM3 %>% filter(CHROM3 == chr3)
    chr3_perm <- perm_test(chr3_LOHcounts_ID, 
                          cat_var = "Tx_name", 
                          cat_names = with(Tx_combo[tx, ], c(Tx_1, Tx_2)), 
                          response_var = "n_LOH", 
                          alpha = 0.05,
                          n_perms = 10000, alt_hyp = "two-tailed")
    all_chr3_perm[c_i, 3:6] <- unlist(chr3_perm[1:4])
    print(paste(all_chr3_perm[1, "Combo"], chr3, sep = " "))
  }
  all_CHROM3_LOHrate_perm <- rbind(all_CHROM3_LOHrate_perm, all_chr3_perm)
}

i_false_sig <- all_CHROM3_LOHrate_perm$pVal >= 0.05 & all_CHROM3_LOHrate_perm$rejectNull == 1
all_CHROM3_LOHrate_perm$rejectNull[i_false_sig] <- 0

all_CHROM3_LOHrate_perm <- BHcorrection(all_CHROM3_LOHrate_perm)

all_CHROM3_LOHrate_perm$sig_lab_p <- ifelse(all_CHROM3_LOHrate_perm$rejectNull == 1, "*", "")
all_CHROM3_LOHrate_perm$sig_lab_BH <- ifelse(all_CHROM3_LOHrate_perm$rejectNull_BH == 1, "**", "")
all_CHROM3_LOHrate_perm$sig_lab_both <- paste0(all_CHROM3_LOHrate_perm$sig_lab_p, 
                                              ifelse(all_CHROM3_LOHrate_perm$rejectNull_BH == 1, ", ", ""),
                                              all_CHROM3_LOHrate_perm$sig_lab_BH)

combo_cols <- colsplit(all_CHROM3_LOHrate_perm$Combo, "-", c("Tx1", "Tx2"))
all_CHROM3_LOHrate_perm <- cbind(all_CHROM3_LOHrate_perm, combo_cols)

all_CHROM3_LOHrate_perm$Tx1 <- factor(all_CHROM3_LOHrate_perm$Tx1, levels = c("WT", "Cas9", "Drive"))
all_CHROM3_LOHrate_perm$Tx2 <- factor(all_CHROM3_LOHrate_perm$Tx2, levels = c("WT", "Cas9", "Drive"))
all_CHROM3_LOHrate_perm$Combo <- factor(all_CHROM3_LOHrate_perm$Combo)
all_CHROM3_LOHrate_perm$CHROM <- factor(substr(all_CHROM3_LOHrate_perm$CHROM, 1, 2))
all_CHROM3_LOHrate_perm$CHROM3_pos <- factor(substr(all_CHROM3_LOHrate_perm$CHROM, 4, 4))
sig_CHROM3rate <- all_CHROM3_LOHrate_perm %>% filter(rejectNull == 1)

y_loc <- max(LOHcounts_ID_CHROM3_mean$mean_rate + LOHcounts_ID_CHROM3_mean$CI_rate)
as.n <- function(x) {
  y <- as.numeric(x)
  return(y)
}

LOHrate_CHROM3_line <- LOHcounts_ID_CHROM3_mean %>% 
  ggplot() + 
  facet_wrap(~CHROM, ncol = 8, strip.position = "bottom", scales = "free_x") +
  geom_errorbar(aes(x = CHROM3_pos, ymin = mean_rate, ymax = mean_rate, color = Tx_name), 
                width = 0.6, position = position_dodge(width = 0.85), size = 0.75) +
  geom_errorbar(aes(x = CHROM3_pos, ymin = ifelse(mean_rate - CI_rate < 0, 0, mean_rate - CI_rate), 
                    ymax = mean_rate + CI_rate, color = Tx_name),
                width = 0, position = position_dodge(width = 0.85), size = 0.75) +
  # geom_bracket(data = sig_CHROM3rate,
  #              y.position = y_loc*1.05 + (as.numeric(sig_CHROM3rate$Combo)-1)*y_loc*0.07,
  #              step.increase = 0, # step.group.by = sig_CHROM3rate$CHROM,
  #              vjust = 0, label.size = 5,
  #              aes(xmin = as.n(CHROM3_pos) + 0.333*0.85*(as.n(Tx1)-2),
  #                  xmax = as.n(CHROM3_pos) + 0.333*0.85*(as.n(Tx2)-2),
  #                  label = sig_lab_both)) +
  geom_text(aes(x = CHROM3_pos, y = -0.00005, label = total_LOH, color = Tx_name),
            position = position_dodge(width = 0.85), show.legend = FALSE, size = 3.5) +
  scale_color_manual(values = txPal, name = "Drive Type") +
  # ylim(c(0, NA)) +
  ylab("LOH event rate (/genome/generation)") + 
  xlab("Chromosome") +
  theme(axis.text.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.background = element_rect(color = "grey80"),
        text = element_text(size = 14))

LOHrate_CHROM3_line

ggsave(file.path(outIntDir, "LOHrate_CHROM3_line_2021_12.png"), 
       plot = LOHrate_CHROM3_line,
       device = "png",
       width = 16, height = 9, 
       units = "in",
       dpi = 600)

# LOHcounts_ID_CHROM3_mean %>% 
#   ggplot(aes(x = chr3, y  = mean_rate, color = Tx_name)) + 
#   geom_errorbar(aes(ymin = ifelse(mean_rate - se_rate < 0, 0, mean_rate - se_rate), 
#                     ymax = mean_rate + se_rate, color = Tx_name)) +
#   geom_point() +
#   # geom_line(aes(group = Tx_name)) +
#   scale_color_manual(values = txPal) + 
#   scale_fill_manual(values = txPal) + 
#   facet_wrap(~CHROM, ncol = 4, scales = "free_x")


LOHcounts_ID_CHROM3_mean %>% ggplot() + geom_histogram(aes(x = log(mean_rate)), bins = 5)

LOHcounts_ID_CHROM3_mean$ord  <- cut(LOHcounts_ID_CHROM3_mean$mean_rate, 
             breaks = nrow(LOHcounts_ID_CHROM3_mean), 
             include.lowest=TRUE,
             labels = str_pad(1:nrow(LOHcounts_ID_CHROM3_mean), width = 3, pad = "0"), ordered=TRUE)

LOHcounts_ID_CHROM3_mean %>% group_by(Tx_name) %>% summarize(v = var(mean_rate))
LOHcounts_ID_CHROM3_mean %>% group_by(chr3) %>% summarize(v = var(mean_rate)) %>% plot()

chr3_anova <- aov(mean_rate ~ chr3, data = LOHcounts_ID_CHROM3_mean)
summary(chr3_anova)

chr3_Tx_anova <- aov(mean_rate ~ Tx_name, data = LOHcounts_ID_CHROM3_mean)
summary(chr3_Tx_anova)

## fit ordered logit model and store results 'm'
m <- kruskal.test(mean_rate ~ Tx_name, data = LOHcounts_ID_CHROM3_mean)
m2 <- kruskal.test(mean_rate ~ chr3, data = LOHcounts_ID_CHROM3_mean)

m3 <- pairwise.wilcox.test(LOHcounts_ID_CHROM3_mean$mean_rate, LOHcounts_ID_CHROM3_mean$Tx_name)
m4 <- pairwise.wilcox.test(LOHcounts_ID_CHROM3_mean$mean_rate, 
                           LOHcounts_ID_CHROM3_mean$chr3, p.adjust.method = "none")

## view a summary of the model
summary(m)

###############################################################################

all_LOHrates_SW50k <- POSi_data_wHet %>% 
  LOHrateSW(by_GT = F, win = 50000, chrom_df = chrom_bound_BY)
mean_LOHrates50k <- all_LOHrates_SW50k %>% 
  rateMeans_xTx(per_gen = T)
mean_LOHrates50k <- mean_LOHrates50k %>% 
  filter(n_clones >= 20)
mean_LOHrates50k <- ChromosomeCoordinates(mean_LOHrates50k)

all_LOHrates_SW50k_I <- POSi_data_wHet %>% 
  filter(!isTerm) %>% 
  LOHrateSW(by_GT = F, win = 50000, chrom_df = chrom_bound_BY)
mean_LOHrates50k_I <- all_LOHrates_SW50k_I %>% 
  rateMeans_xTx(per_gen = T)
mean_LOHrates50k_I <- mean_LOHrates50k_I %>% 
  filter(n_clones >= 20)
mean_LOHrates50k_I <- ChromosomeCoordinates(mean_LOHrates50k_I)

all_LOHrates_SW50k_T <- POSi_data_wHet %>% 
  filter(isTerm)) %>% 
  LOHrateSW(by_GT = F, win = 50000, chrom_df = chrom_bound_BY)
mean_LOHrates50k_T <- all_LOHrates_SW50k_T %>%
  rateMeans_xTx(per_gen = T) %>% 
  filter(n_sites > 0 & !is.na(rate))
mean_LOHrates50k_T <- mean_LOHrates50k_T %>% 
  filter(n_clones >= 20)
mean_LOHrates50k_T <- ChromosomeCoordinates(mean_LOHrates50k_T)

smooth_sw_in <- mean_LOHrates50k # %>% filter(n_sites >= 15000)
smooth_sw_in <- mean_LOHrates50k_I
smooth_sw_in$per_gen_rate <- smooth_sw_in$rate
smooth_sw_in$per_gen_se <- smooth_sw_in$std_err

smooth_sw_in$end_POS <- ConvertPosIndicies(smooth_sw_in, pos_col = "end_POSi", index_out = "POS")

smooth_sw_in$plot_POS <- (smooth_sw_in$POS + smooth_sw_in$end_POS)/2

lims <- smooth_sw_in %>% 
  group_by(CHROM) %>% 
  summarise(xmin=min(POS), xmax=max(end_POS))

lims$breaks <- floor((lims$xmax - lims$xmin)/50000)

ylims <- c(0, NA)

SW_LOH_plot <- ggplot() + 
  geom_rect(data = chrom_arms,
            aes(xmin = 0, xmax = (abs(arm_1_cent) + arm_2_cent)/1000,
                ymin = -Inf, ymax = Inf),
            alpha = 0.2, fill = "white", color = "grey70") +
  geom_vline(data = chrom_arms,
             aes(xintercept = abs(arm_1_cent)/1000), size = 0.35, color = "blue4", alpha = 0.5) +
  # geom_point(data = centrom_df[iCh, ], aes(x = POS, y = 0), color = "blue4", shape = 17, size = 1.5) +
  geom_ribbon(data = smooth_sw_in, 
              aes(x = plot_POS/1000, 
                  y = per_gen_rate,
                  xmin = plot_POS/1000, xmax = plot_POS/1000, 
                  ymin = ifelse((per_gen_rate - per_gen_se) < 0, 0, (per_gen_rate - per_gen_se)), 
                  ymax = (per_gen_rate + per_gen_se),
                  fill = Tx_name),
              alpha = 0.2) +
  geom_line(data = smooth_sw_in,
            aes(x = plot_POS/1000, y = per_gen_rate, color = Tx_name),
            size = 0.35)  +
  geom_point(data = smooth_sw_in,
             aes(x = plot_POS/1000, y = per_gen_rate, color = Tx_name),
             size = 0.2)  +
  # geom_blank(data = subset(smooth_sw_in, CHROM == "04"), 
  #            aes(x = plot_POS, y = rate/n_gens, color = Tx_name)) +
  scale_color_manual(values = txPal, name = "Drive type") +
  scale_fill_manual(values = txPal, name = "Drive type") +
  scale_x_continuous(expand = expansion(0, 0)) +
  scale_y_continuous(expand = expansion(0, 5E-6), 
                     labels = function(x) format(x, scientific = TRUE),
                     limits = ylims
                     # breaks = c(0, 1E-4, 2E-4, 3E-4, 4E-4, 5E-4)
  ) +
  theme(panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(color = "grey40", size = 0.1),
        legend.position = "bottom") +
  xlab("Chromosome position (kbp)") + ylab("LOH rate (/bp/generation)") +
  facet_wrap(~CHROM, ncol = 4, strip.position = "top", scales = "free_x") 

SW_LOH_plot

ggsave(file.path(outIntDir, "bCHROM_LOHrate_plot_2021_11.png"), 
       plot = SW_LOH_plot,
       device = "png",
       width = 16, height = 9, 
       units = "in",
       dpi = 600)

###############################################################################
# LOH rate vs chromsome size

LOHcounts_ID_CHROM_mean$chrom_length <- rep(chrom_lengths_BY, 3)

LOH_chrom_reg <- lm(mean_LOH ~ chrom_length, data = LOHcounts_ID_CHROM_mean)
summary(LOH_chrom_reg)
LOH_chrom_reg$coefficients


LOHcounts_ID_CHROM_mean %>% ggplot() +
  stat_smooth(data = LOH_chrom_reg$model, 
              aes_string(x = names(LOH_chrom_reg$model)[2], y = names(LOH_chrom_reg$model)[1]), 
              method = "lm", col = "grey20", size = 0.5) +
  geom_point(aes(x = chrom_length, y = mean_LOH, color = Tx_name)) +
  # geom_abline(aes(intercept = LOH_chrom_reg$coefficients[1], slope = LOH_chrom_reg$coefficients[2])) +
  scale_color_manual(values = txPal) +
  ylim(0, NA) + xlim(0,NA) +
  geom_text(check_overlap = T, hjust = 0, aes(x = 1.1E6, y = 0.15), 
            label = paste0("y = ", formatC(LOH_chrom_reg$coef[[2]], 2, format = "E"), 
                           "x + ", signif(LOH_chrom_reg$coef[[1]], 3),
                          "\nAdj R2 = ", signif(summary(LOH_chrom_reg)$adj.r.squared, 3),
                     # "\nIntercept =",signif(LOH_chrom_reg$coef[[1]], 3),
                     # "\nSlope =", formatC(LOH_chrom_reg$coef[[2]], 2, format = "E"),
                     "\np = ",formatC(summary(LOH_chrom_reg)$coef[2,4], 2, format = "E")))
  # facet_wrap(~Tx_name, ncol = 1)

ggplotReg <- function (fit) {
  
  require(ggplot2)
  
  ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + 
    geom_point() +
    stat_smooth(method = "lm", col = "red") +
    labs(title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
                       "Intercept =",signif(fit$coef[[1]],5 ),
                       " Slope =",signif(fit$coef[[2]], 5),
                       " P =",signif(summary(fit)$coef[2,4], 5)))
}

ggplotReg(LOH_chrom_reg)


#!!!###########################################################################
# Rate of breakpoints per clone per position window within chromosomes

# POSi_data_in$est_mid_crct <- ConvertPosIndicies(POSi_data_in, pos_col = "est_mid_POS", index_out = "POSi")

POSi_data_in <- POSi_data_in %>% mutate(bp = 1)
LOH_mid_SW_coarse <- POSi_data_in %>% 
  filter(!isTerm) %>%
  # filter(!cross_h_term) %>% 
  SliderCalc(data_col = "bp", index_col = "est_mid", 
             window_size = 50000, slide_interval = 1000,
             summary_stat = sum, factor_col = "Tx_name", chrom_win = T)
LOH_mid_SW_coarse$Tx_name <- factor(LOH_mid_SW_coarse$Tx_name, levels = Tx_name_levels)

LOH_mid_SW_fine <- POSi_data_in %>% 
  filter(!isTerm) %>%
  # filter(!cross_h_term) %>% 
  SliderCalc(data_col = "bp", index_col = "est_mid", 
             window_size = 100000, slide_interval = 100,
             summary_stat = sum, factor_col = "Tx_name", chrom_win = T)
LOH_mid_SW_fine$Tx_name <- factor(LOH_mid_SW_fine$Tx_name, levels = Tx_name_levels)

LOH_ID_SW <- POSi_data_in %>% 
  filter(!isTerm) %>%
  # filter(!cross_h_term) %>% 
  SliderCalc(data_col = "bp", index_col = "est_mid", 
             window_size = 50000, slide_interval = 50000,
             summary_stat = sum, factor_col = "ID", chrom_win = T)

LOH_ID_SW$Tx <- factor(substr(LOH_ID_SW$ID, 1, 1))
LOH_ID_SW$Tx_name <- Recode_Tx(LOH_ID_SW$Tx)
LOH_ID_SW <- LOH_ID_SW %>% select(-Tx)

fract_chrom_fun <- function(df, reps = 10000, alpha = 0.05) {
  # df <- LOH_ID_SW
  t_lo <- floor(reps * (alpha/2))
  t_up <- ceiling(reps * (1 - (alpha/2)))
  CI_df <- data.frame(NULL)
  for(tx in Tx_name_levels) {
    # tx <- "WT"
    n_id <- n_clones_xTx %>% filter(Tx_name == tx) %>% pull(n)
    tx_CI_df <- data.frame(NULL)
    for(chr in chrom_IDs$CHROM) {
      # chr <- "02"
      sub_df <- df %>% filter(Tx_name == tx, CHROM == chr)
      s_POSi <- sub_df %>% distinct(start) %>%
        arrange(start) %>% pull(start)
      n_s <- length(s_POSi)
      fract_m <- matrix(nrow = n_s, ncol = reps)
      for(rep in 1:reps) {
        # rep <- 8
        sample_df <- sub_df %>% group_by(start) %>% sample_n(size = n_id, replace = T)
        n_loh <- sum(sample_df$bp)
        n_loh <- ifelse(n_loh == 0, 1, n_loh)
        rep_v <- sample_df %>% group_by(start) %>% summarize(f = sum(bp)/n_loh) %>% pull(f)
        fract_m[, rep] <- rep_v
        # colnames(rep_df)[2] <- paste0("r_", rep)
        # fract_df <- cbind(fract_df, rep_df)
        
      }
      fract_m_sort <- t(apply(fract_m, 1, sort))
      chrom_CI_df <- data.frame(Tx_name = tx, CHROM = chr, start = s_POSi, 
                          CI_lo = fract_m_sort[, t_lo], CI_up = fract_m_sort[, t_up])
      tx_CI_df <- rbind(tx_CI_df, chrom_CI_df)
    }
    CI_df <- rbind(CI_df, tx_CI_df)
  }
  CI_df$Tx_name <- factor(CI_df$Tx_name, levels = Tx_name_levels)
  CI_df <- CI_df %>% arrange(Tx_name, start)
  return(CI_df)
}

LOH_SW_CIs <- fract_chrom_fun(LOH_ID_SW)
iLOH_SW_CIs <- LOH_SW_CIs

LOH_mid_SW <- LOH_mid_SW_coarse
# colnames(LOH_SW_CIs)[grep("start", colnames(LOH_SW_CIs))] <- "start"
LOH_CHROM_counts <- POSi_data_in %>% filter(GT != "0/1", !isTerm) %>% 
  group_by(Tx_name, CHROM) %>% summarize(tot_LOH = n())
for(tx in Tx_name_levels) {
  # tx = "WT"
  # x = sub_counts$CHROM[1]
  i_SW <- LOH_mid_SW$Tx_name == tx
  sub_df <- LOH_mid_SW %>% filter(Tx_name == tx)
  sub_counts <- LOH_CHROM_counts %>% filter(Tx_name == tx)
  fract_list <- sapply(sub_counts$CHROM, function(x) 
    sub_df$bp[sub_df$CHROM == x]/sub_counts$tot_LOH[sub_counts$CHROM == x])
  LOH_mid_SW$chrom_fraction[i_SW] <- unlist(fract_list)
}


# LOH_mid_SW <- LOH_mid_SW %>% group_by(Tx_name, CHROM) %>% mutate(chrom_fraction = bp/sum(bp))

# LOH_mid_SW <- merge(LOH_mid_SW, LOH_SW_CIs, by = c("Tx_name", "CHROM", "start")) %>% arrange(Tx_name, start)

# Per clone rate. Each Tx has a different # of clones
# LOH_mid_SW$clone_rate <- 0
# LOH_mid_SW$clone_rate[LOH_mid_SW$Tx_name == "WT"] <- 
#   LOH_mid_SW$bp[LOH_mid_SW$Tx_name == "WT"]/n_evo_xTx$n[n_evo_xTx$Tx_name == "WT"]
# LOH_mid_SW$clone_rate[LOH_mid_SW$Tx_name == "Cas9"] <- 
#   LOH_mid_SW$bp[LOH_mid_SW$Tx_name == "Cas9"]/n_evo_xTx$n[n_evo_xTx$Tx_name == "Cas9"]
# LOH_mid_SW$clone_rate[LOH_mid_SW$Tx_name == "Drive"] <- 
#   LOH_mid_SW$bp[LOH_mid_SW$Tx_name == "Drive"]/n_evo_xTx$n[n_evo_xTx$Tx_name == "Drive"]

LOH_mid_SW <- LOH_mid_SW %>% mutate(start_POSi = round(start), end_POSi = round(end)) %>% select(!start:end)
s_POS <- LOH_mid_SW %>% ConvertPosIndicies(pos_col = "start_POSi", index_out = "POS", add_chroms = T)
names(s_POS)[1] <- "start_POS"
e_POS <- LOH_mid_SW %>% ConvertPosIndicies(pos_col = "end_POSi", index_out = "POS")
LOH_mid_SW <- data.frame(LOH_mid_SW, s_POS[c(2,1)], end_POS = e_POS)
colnames(LOH_mid_SW)[1] <- c("Tx_name")
LOH_mid_SW$Tx_name <- factor(LOH_mid_SW$Tx_name, levels = c("WT", "Cas9", "Drive"))
LOH_mid_SW <- LOH_mid_SW %>% 
  mutate(mid_POS = start_POS + 50000/2, mid_POSi = (start_POSi + end_POSi)/2,
         mid_POS_kb = (start_POS + 50000/2)/1000, mid_POSi_kb = (start_POSi + end_POSi)/2000)

LOH_mid_SW$rom_CHROM <- factor(chrom_IDs$rom_CHROM[as.numeric(LOH_mid_SW$CHROM)], levels = roman_chr)

breakpoints_SW_plot <- LOH_mid_SW %>% 
  # filter(CHROM == "03") %>%
  # filter(Tx_name == "Drive") %>%
  ggplot(aes(x = mid_POS_kb, y = chrom_fraction)) + 
  geom_vline(data = centrom_df, aes(xintercept = POS/1000), 
             linetype = 2,color = "blue4", size = 0.25) +
  # geom_ribbon(aes(ymin = CI_lo, ymax = CI_up, fill = Tx_name, color = Tx_name), 
  #             size = 0.05, alpha = 0.1) +
  geom_line(aes(color = Tx_name)) +
  # geom_line(data = genome_scores_50kb, 
  #           aes(x = mid_POS/1000, y = -(score - min_score)/((max(score) - min_score)*3))) +
  # geom_line(aes(x = mid_POS_kb, y = log(rate, 10), color = Tx_name)) +
  # geom_point(aes(x = mid_POS_kb, y = rate, color = Tx_name)) +
  # xlim(0, max(chrom_lengths_BY)) +
  xlab("Window Position (kb)") + ylab("Breakpoints in window / Breakpoints in chromosome") +
  xlim(0,NA) +
  # scale_y_continuous(limits = c(-1/3, 0.4), 
  #                    breaks = c(-1/3, -2/9, -1/9, seq(0, 0.4, 0.1)), 
  #                    labels = c(210, 140, 70, seq(0, 0.4, 0.1))) +
  scale_color_manual(values = txPal, name = "Strain") + 
  scale_fill_manual(values = txPal, name = "Strain") + 
  facet_wrap(~rom_CHROM, ncol = 4, scales = "free_x") +
  theme(legend.text = element_text(size = 13),
        legend.title = element_text(size = 14),
        strip.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 16),
        # text = element_text(size = 18), 
        plot.margin = margin(20, 20, 20, 20),
        panel.grid.minor = element_blank(),
        legend.position = c(0.97, 0.94),
        legend.background = element_rect(fill = "white", color = "white"))
  
breakpoints_SW_plot

ggsave(file.path(outIntDir, "iLOH_bp_SW50x1_plot.png"), 
       plot = breakpoints_SW_plot,
       device = "png",
       width = 16, height = 9, 
       units = "in",
       dpi = 600)

POSi_data_in %>% 
  ggplot() + geom_histogram(aes(x = est_mid, fill = Tx_name), binwidth = 10000) +
  # geom_point(aes(x = mid_POS, y = value, color = Tx_name)) +
  # xlim(0, max(chrom_lengths_BY)) +
  xlab("Window Position (kbp)") + ylab("Number of Breakpoints") +
  xlim(0, NA) +
  scale_fill_manual(values = txPal, name = "Drive Type") + 
  facet_wrap(~CHROM, ncol = 4, scales = "free_x") +
  theme(text = element_text(size = 14), legend.position = "bottom")

DrivevsWT_breaks <- LOH_mid_SW %>% 
  select(Tx_name, CHROM, rate, mid_POS_kb, mid_POSi_kb) %>%
  filter(Tx_name != "Cas9") %>% 
  pivot_wider(names_from = Tx_name, values_from = rate) %>%
  mutate(rate_diff = (Drive - WT)/(Drive + WT + 0.00001))

genome_scores_50kb <- genome_scores_50kb %>%
  mutate(mid_POS_kb = mid_POS/1000, norm_score = (score - min_score)/(max(score) - min_score))

DrivevsWT_POS_plot <- DrivevsWT_breaks %>% 
  # filter(CHROM == "05") %>%
  ggplot() + 
  geom_line(aes(x = mid_POS_kb, y = rate_diff)) +
  geom_line(data = genome_scores_50kb,
            aes(x = mid_POS_kb, y = norm_score), color = "red3") +
  # geom_line(aes(x = mid_POS_kb, y = log(rate, 10), color = Tx_name)) +
  # geom_point(aes(x = mid_POS_kb, y = rate, color = Tx_name)) +
  # xlim(0, max(chrom_lengths_BY)) +
  xlab("Window Position (kbp)") + ylab("Number of Breakpoints / Clone") +
  xlim(0,NA) +
  # scale_y_continuous(limits = c(-1/3, 0.3), 
  #                    breaks = c(-1/3, -2/9, -1/9, seq(0, 0.3, 0.1)), 
  #                    labels = c(210, 140, 70, seq(0, 0.3, 0.1))) +
  # scale_color_manual(values = txPal, name = "Drive Type") + 
  facet_wrap(~CHROM, ncol = 4, scales = "free_x") +
  theme(text = element_text(size = 14), legend.position = "bottom")

DrivevsWT_POS_plot

DrivevsWT_breaks_scores <- cbind(DrivevsWT_breaks, norm_score = genome_scores_50kb$norm_score)

DrivevsWT_breaks_scores %>% 
  # filter(!(WT == 0 & Drive == 0)) %>%
  ggplot() + 
  geom_point(aes(x = rate_diff, y = norm_score)) +
  xlab("Breakpoint rate bias") + ylab("Normailzed gRNA similarity score") +
  ylim(0,1) +
  scale_x_continuous(breaks = c(-1, -0.5, 0, 0.5, 1), 
                     labels = c("WT only", "WT higher", "0", "Drive higher", "Drive only")) + 
  theme(panel.grid.minor.x = element_blank())


  
    