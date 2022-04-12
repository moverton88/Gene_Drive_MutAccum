# Get set of SNP mutations from SNPs_finalGT


##########

# denovo_SNPs <- denovo_SNPs %>% filter(!Cut, GQ >= 50, !existing_SNP)

dn_GT_vals <- denovo_SNPs %>% site_genotype_stats()

# Get sitewise genotype stats for each founder group
dn_Line_GT_vals <- data.frame(NULL)
for(l in levels(denovo_SNPs$Line)) {
  # l <- "N_C"
  line_SNPs <- denovo_SNPs %>% filter(Line == l)
  line_SNPs$Line <- droplevels(line_SNPs$Line)
  line_GT_vals <- line_SNPs %>% site_genotype_stats()
  line_GT_vals$Line <- l
  dn_Line_GT_vals <- rbind(dn_Line_GT_vals, line_GT_vals)
}

# Get sites for which there is a single polymorphism among all clones of a founder group
evo_single_mut_Line <- dn_Line_GT_vals %>% 
  filter(nRef_evo >= 7 & ((nHet_evo == 1 & nAlt_evo == 0) | (nAlt_evo == 1 & nHet_evo == 0))) %>% 
  select(Line, POSi, nRef_anc, nRef_evo, nHet_evo)

# Get positions were the above is true
evo_single_Line_POSi <- evo_single_mut_Line %>%
  distinct(POSi) %>% pull(POSi)

# Find sites for which there is a single polymorphism among all end-point clones and get postions
evo_single_mut_POSi <- dn_GT_vals %>% 
  filter(POSi %in% evo_single_Line_POSi & 
           ((nHet_evo == 1 & nAlt_evo == 0) | (nHet_evo == 0 & nAlt_evo == 1))) %>%
  pull(POSi)

# Sites for which all founders are homozygous reference
anc_hom_POSi <- dn_GT_vals %>% filter(fRef_anc == 1) %>% pull(POSi) %>% unique()

# evo_hom_POSi <- dn_GT_vals %>% filter(nRef_evo >= 30) %>% pull(POSi) %>% unique()

# Sites for which only a single end-point clone is heterozygous
evo_single_het_POSi <- dn_GT_vals %>% 
  filter(nHet_evo == 1 & nAlt_evo == 0) %>% pull(POSi) %>% unique()

# Sites for which only a single end-point clone is homozygous alt
evo_single_alt_POSi <- dn_GT_vals %>% 
  filter(nHet_evo == 0 & nAlt_evo == 1) %>% pull(POSi) %>% unique()

# Sites for which two end-point clones are heterozygous
evo_double_het_POSi <- dn_GT_vals %>% 
  filter(nHet_evo == 2 & nAlt_evo == 0) %>% pull(POSi) %>% unique()

# Unique heterozygous variant sites
anc_hom_mut_het_POSi <- intersect(anc_hom_POSi, evo_single_het_POSi)
# evo_hom_mut_het_POSi <- intersect(evo_hom_POSi, evo_single_het_POSi)
# mut_het_POSi <- union(anc_hom_mut_het_POSi, evo_hom_mut_het_POSi)

# Unique homozygous variant sites
anc_hom_mut_alt_POSi <- intersect(anc_hom_POSi, evo_single_alt_POSi)
# evo_hom_mut_alt_POSi <- intersect(evo_hom_POSi, evo_single_alt_POSi)
# mut_alt_POSi <- union(anc_hom_mut_alt_POSi, evo_hom_mut_alt_POSi)

# mut_all_POSi_2 <- union(mut_het_POSi, mut_alt_POSi)

# Unique variant sites
mut_all_POSi <- union(union(anc_hom_mut_het_POSi, anc_hom_mut_alt_POSi), evo_single_mut_POSi)

# SNPs_merge_finalGT %>% mutate(single_het = POSi %in% evo_single_het_POSi) %>%
#   filter(!Cut, GQ >= 50, !(GT != "0/0" & Alt_DP_final < 3), Sum_DP_final < 300) %>% 
#   ggplot() + geom_histogram(aes(x = Sum_DP_final)) + facet_wrap(~single_het, ncol = 1, scales = "free_y", as.table = F) +
#   ggtitle("Read depth for sites with single het end-point clones")


# final_mut_het <- denovo_SNPs %>% 
#   filter(POSi %in% mut_all_POSi, Rep != "00", GT == "0/1") %>% arrange(ID, POSi)
# 
# final_mut_alt <- denovo_SNPs %>% 
#   filter(POSi %in% mut_all_POSi_2, Rep != "00", GT == "1/1") %>% arrange(ID, POSi)
# 
# final_mut_set <- rbind(final_mut_het, final_mut_alt) %>% arrange(ID, POSi)


final_mut_set <- denovo_SNPs %>% 
  filter(POSi %in% mut_all_POSi, Rep != "00", GT %in% c("0/1", "1/1")) %>% arrange(ID, POSi)

evo_IDs <- SNPs_merge_finalGT %>% distinct(ID, .keep_all = T) %>% filter(Rep != "00") %>% pull(ID)
final_mut_set$ID <- factor(final_mut_set$ID, levels = evo_IDs)

final_mut_set <- final_mut_set %>% 
  group_by(ID) %>%
  mutate(diff = POSi - lag(POSi),
         lead_diff = lead(POSi) - POSi)
final_mut_set$diff[is.na(final_mut_set$diff)] <- 12.2E7
final_mut_set$lead_diff[is.na(final_mut_set$lead_diff)] <- 12.2E7

# final_mut_set %>% filter(diff < 100 | lead_diff < 100) %>% select(!QUAL_BYcall:Tx_name) %>% View()

# All clustered SNPs are in Ty LTRs and other types of repeats
final_mut_set_merge <- final_mut_set %>% filter(diff > 100 & lead_diff > 100)
final_mut_set_merge$CHROM <- factor(final_mut_set_merge$CHROM, levels = str_pad(1:16, width = 2, pad = "0"))
n_complex <- nrow(final_mut_set) - nrow(final_mut_set_merge)

# Since triploid founder groups have similar SNM rates, merge with main set
###############################################################################
final_mut_set_merge <- rbind(final_mut_set_merge, triploid_mut_set_merge)
final_mut_set_merge <- final_mut_set_merge %>% arrange(Tx_name, ID, POSi) %>% select(!(GQ:existing_SNP))
final_mut_set_merge$Tx_ID <- Recode_Tx_ID(final_mut_set_merge$Tx)

all_evo_IDs <- c(evo_IDs, as.character(triploid_evo_IDs)) %>% sort()
final_mut_set_merge$ID <- factor(final_mut_set_merge$ID, levels = all_evo_IDs)

###############################################################################
# Mutation rate
# n_clones <- SNPs_merge_finalGT %>% filter(Rep != "00") %>% pull(ID) %>% unique() %>% length()
# n_clones_xTx <- SNPs_merge_finalGT %>% 
#   filter(Rep != "00") %>%
#   distinct(ID, .keep_all = T) %>% count(Tx)

n_SNMs <- nrow(final_mut_set_merge)

# Number of diploid positions
# n_sites <- (mean_cover - n_dropped_sites) * 2
n_sites <- mean_cover * 2
n_clones <- length(levels(final_mut_set_merge$ID))
SNM_rate <- n_SNMs / (n_clones*n_gens*n_sites)

SNM_rate

###############################################################################
# Compare distribution of SNMs among clones
SNMs_final_counts <- final_mut_set_merge %>% count(ID) 

evo_zeros <- data.frame(ID = all_evo_IDs[!all_evo_IDs %in% SNMs_final_counts$ID], n = 0)

SNMs_final_counts <- rbind(evo_zeros, SNMs_final_counts) %>% arrange(ID)

SNMs_final_counts <- SNMs_final_counts %>% CategoriesFromID()


# SNMs_final_counts %>% group_by(Line %in% noAncestor) %>% count(n)
mean_SNMrate <- SNMs_final_counts %>% 
  summarise(mean_SNM = mean(n), 
            sd_SNM = sd(n), n_clones = n()) %>% 
  mutate(se_SNM = sd_SNM/sqrt(n_clones)) %>% 
  mutate(mean_rate = mean_SNM / n_gens / n_sites, 
            se_rate = se_SNM / n_gens / n_sites)

tx_SNMrate <- SNMs_final_counts %>% 
  group_by(Tx_ID) %>% 
  summarise(mean_SNM = mean(n), 
            sd_SNM = sd(n), n_clones = n()) %>% 
  mutate(se_SNM = sd_SNM/sqrt(n_clones)) %>% 
  mutate(mean_rate = mean_SNM / n_gens / n_sites, 
         se_rate = se_SNM / n_gens / n_sites)

# Output table as .png
png(paste0(outIntDir, "SNM_rate_table.png"), height = 50*nrow(tx_SNMrate), width = 100*ncol(tx_SNMrate))
grid.table(mean_SNMrate, theme = ttheme_minimal(), rows = NULL)
dev.off()

# SNMs_final_counts %>%
#   ggplot() + geom_histogram(aes(x = n, y = ..density.., fill = Tx_name), 
#                             binwidth = 1,
#                             position = position_dodge(preserve = "single", width = 0.98)) +
#   scale_x_continuous(breaks = seq(0, 6, 1)) +
#   scale_fill_manual(values = txPal)

SNMs_final_counts_table <- SNMs_final_counts %>% 
  group_by(Tx_name) %>% count(n, name = "n_muts") %>% 
  arrange(Tx_name)

colnames(SNMs_final_counts)[which(colnames(SNMs_final_counts) == "n_SNMs")] <- "n"

SNMs_final_counts_stats <- SNMs_final_counts %>% 
  group_by(Tx_ID) %>% 
  summarize(n_clones = n(),
            mean_n_PM = mean(n), se_n = se(n), 
            meanCI95lo = (mean(n) - se(n)*1.96), 
            meanCI95hi = (mean(n) + se(n)*1.96),
            bprate = mean(n)/n_gens/(mean_cover*2), 
            bprateCI95lo = (mean(n) - se(n)*1.96)/n_gens/(mean_cover*2), 
            bprateCI95hi = (mean(n) + se(n)*1.96)/n_gens/(mean_cover*2))

# the count function does not include 0 counts, so we must do this manually

SNMs_final_counts_table <- fill_zeros(SNMs_final_counts_table)

SNMs_final_counts_table <- SNMs_final_counts_table %>% group_by(Tx_name) %>%
  mutate(density = n_muts/sum(n_muts)) %>% ungroup()

SNMs_final_counts_table$Tx_name <- factor(SNMs_final_counts_table$Tx_name, 
                                          levels = c("WT", "Cas9", "Drive"))

SNMs_means <- SNMs_final_counts %>% group_by(Tx_name) %>% summarize(m = mean(n), s = sd(n))
# SNMs_means$m[1] <- SNMs_means$m[1] + 0.01
# SNMs_means$m[2] <- SNMs_means$m[2] - 0.01

SNMs_counts_table_wide <- SNMs_final_counts_table %>% 
  select(-c(density)) %>% 
  pivot_wider(names_from = Tx_name, values_from = n_muts) %>% 
  as.data.frame() %>% arrange(n_SNMs)

SNMs_counts_table_wide[is.na(SNMs_counts_table_wide)] <- 0
N_counts <- SNMs_counts_table_wide[, "WT"]
H_counts <- SNMs_counts_table_wide[, "Cas9"]
F_counts <- SNMs_counts_table_wide[, "Drive"]
i_WT <- SNMs_counts_table_wide %>% filter(WT != 0) %>% nrow()

i_pos <- !SNMs_counts_table_wide[,-1] == 0

mut_WT_Cas9_chisq <- chisq.test(cbind(N_counts[i_pos[, 1] | i_pos[, 2]], H_counts[i_pos[, 1] | i_pos[, 2]]))

mut_WT_Drive_chisq <- chisq.test(cbind(N_counts[i_pos[, 1] | i_pos[, 3]], F_counts[i_pos[, 1] | i_pos[, 3]]))

mut_Cas9_Drive_chisq <- chisq.test(cbind(H_counts[i_pos[, 2] | i_pos[, 3]], F_counts[i_pos[, 2] | i_pos[, 3]]))

mut_WT_Cas9_perm <- perm_test(SNMs_final_counts, 
                              cat_var = "Tx_name", cat_names = c("WT", "Cas9"), 
                              response_var = "n", alpha = 0.05)

mut_WT_Drive_perm <- perm_test(SNMs_final_counts, 
                              cat_var = "Tx_name", cat_names = c("WT", "Drive"), 
                              response_var = "n", alpha = 0.05)

mut_Cas9_Drive_perm <- perm_test(SNMs_final_counts, 
                              cat_var = "Tx_name", cat_names = c("Cas9", "Drive"), 
                              response_var = "n", alpha = 0.05)


SNM_counts_label <- paste0("Total SNMs = ", n_SNMs, 
                           "\nPermutation test \n  WT - Cas9 p = ", round(mut_WT_Cas9_perm$p_value, 3),
                           "\n  WT - Drive p = ", round(mut_WT_Drive_perm$p_value, 3),
                           "\n  Cas9 - Drive p = ", round(mut_Cas9_Drive_perm$p_value, 3)
                           )

lab_df <- data.frame(Tx_name = factor("Drive", levels = c("WT", "Cas9", "Drive")), 
                     n_SNMs = 5, n_muts = 20)

# Dotplot of point mutation count !!!!!!!######################################
SNMs_final_counts$Tx_ID <- substr(SNMs_final_counts$Tx_name, 1, 1)
SNMs_final_counts$Tx_ID <- factor(SNMs_final_counts$Tx_ID, levels = c("D", "C", "W"))
SNMs_final_counts <- SNMs_final_counts %>% arrange(Tx_ID, n)
DA_clones <- c("F_A09", "F_C02", "F_D01", "F_F03", "F_F07", "F_G10")

SNMs_final_counts$dot_color <- txPal[as.numeric(SNMs_final_counts$Tx_name)]
SNMs_final_counts$fill_color <- txPal[as.numeric(SNMs_final_counts$Tx_name)]
SNMs_final_counts$fill_color[SNMs_final_counts$ID %in% DA_clones] <- NA
SNMs_final_counts_stats$fill_color <- rev(txPal[as.numeric(SNMs_final_counts_stats$Tx_ID)])

lab_y <- max(SNMs_final_counts$n)
nudge <- c(2, 6, 10)
max_break <- ceiling(lab_y/5)*5

sr <- 1.1
x_diff <- (sr - 1)/50

SNMcount_dotplot_mean <- SNMs_final_counts %>%
  # filter(Tx_name == "Cas9") %>%
  ggplot(aes(x = x_diff)) + 
  # geom_rect(data = SNMs_final_counts_stats, 
  #           aes(xmin = -0.45, xmax = 0.45, 
  #               ymin = meanCI95lo, ymax = meanCI95hi),
  #               fill = SNMs_final_counts_stats$fill_color, 
  #               color = NA,
  #           alpha = 0.3) +
  geom_dotplot(aes(y = n, group = Tx_ID), 
               fill = SNMs_final_counts$fill_color, 
               color = SNMs_final_counts$dot_color,
               stroke = 2,
               binwidth = 1, binaxis = "y", 
               stackdir = "center", binpositions = "all",
               dotsize = 0.5, stackratio = sr) +
  stat_summary(aes(y = n), color = "grey30", fun = "mean", fun.min = "mean", fun.max= "mean",
               size = 0.2, width = 0.9, geom = "crossbar") +
  scale_fill_manual(values = txPal) +
  # scale_color_manual(values = txPal) +
  scale_y_continuous(breaks = seq(0, max_break, 5),
                     name = "Number of point mutations") +
  scale_x_continuous(breaks = 0) +
  theme(axis.ticks.y = element_blank(), 
        # axis.ticks.x = element_line(size = 0.25, color = "grey80"),
        # axis.ticks.length = unit(4, "mm"),
        plot.margin = unit(c(t = 0, r = 5, b = 5, l = 5), "mm"),
        panel.grid.minor.x = element_blank(),
        # panel.grid.major.x = element_blank(),
        # panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        # axis.title.y = element_text(hjust = 0.17, vjust = 1.5),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        strip.text = element_text(size = 30),
        text = element_text(size = 30),
        # axis.text = element_text(size = 12),
        legend.position = "none") +
  # facet_grid(.~Line %in% c("H_F", "H_H"))
  facet_grid(.~Tx_ID)

SNMcount_dotplot_mean

ggsave(file.path(outIntDir, "SNMcount_dotplot_mean_2022_03.png"),
       plot = SNMcount_dotplot_mean,
       device = "png",
       width = 16, height = 9,
       units = "in",
       dpi = 600)


SNM_n_plot <- SNMs_final_counts_table %>%
  ggplot() + 
  geom_vline(data = SNMs_final_counts_stats, aes(xintercept = mean_n, color = Tx_name)) +
  geom_rect(data = SNMs_final_counts_stats, aes(xmin = meanCI95lo, xmax = meanCI95hi, 
                                                ymin = 0, ymax = Inf, fill = Tx_name), alpha = 0.4) +
  # geom_label(data = SNMs_final_counts_stats, aes(x = mean_n, y = 30), 
  #            label = paste0(tx_SNMrate$rate, "/bp/gen"),
  #            hjust = "left", label.size = 0, size = 3.5) +
  geom_col(aes(x = n_SNMs, y = n_muts, fill = Tx_name),
           position = position_dodge(preserve = "single")) +
  # geom_label(data = lab_df, aes(x = n_SNMs, y = n_muts), label = SNM_counts_label,
  #            hjust = "left", label.size = 0, size = 4) +
  scale_x_continuous(breaks = seq(0, 15, 1)) +
  xlab("Number of Point Mutations") + ylab("Number of Clones") +
  scale_fill_manual(values = txPal, name = "Treatment") + 
  scale_color_manual(values = txPal, name = "Treatment") +
  coord_flip() +
  theme(legend.position = "bottom", 
        text = element_text(size = 14),
        strip.text.x = element_blank()) +
  facet_grid(.~Tx_name)

SNM_n_plot

# Plot center-out histogram of point mutation count !!!!!!!####################
n_Tx <- SNMs_final_counts_table %>% 
  group_by(Tx_name) %>% summarize(n_clones = sum(n_muts))
n_muts_lim <- max(SNMs_final_counts_table$n_muts) + 1
y_gap <- 0.05

SNMs_final_counts_table$x_left <- -SNMs_final_counts_table$n_muts/2
SNMs_final_counts_table$x_right <- SNMs_final_counts_table$n_muts/2
SNMs_final_counts_table$y_bot <- SNMs_final_counts_table$n_SNMs - 1/2 + y_gap
SNMs_final_counts_table$y_top <- SNMs_final_counts_table$n_SNMs + 1/2 - y_gap

SNMs_final_counts_table$text_pos <- ifelse(SNMs_final_counts_table$n_muts == 0 |
                                             SNMs_final_counts_table$n_muts > 2, 
                                           0, SNMs_final_counts_table$x_right + 2)

SNMs_final_counts_table$text_color <- ifelse(SNMs_final_counts_table$n_muts > 2, 
                                             "white",
                                           txPal[as.numeric(SNMs_final_counts_table$Tx_name)])

text_backing <- SNMs_final_counts_table %>% 
  filter(text_color != "white") %>% select(Tx_name, n_SNMs, text_pos)
text_backing$x_min <- text_backing$text_pos - 1.75
text_backing$x_max <- text_backing$text_pos + 1.75
text_backing$y_min <- text_backing$n_SNMs - 0.3
text_backing$y_max <- text_backing$n_SNMs + 0.3

text_backing_2 <- SNMs_final_counts_table %>% filter(Tx_name == "Cas9", n_SNMs == 1) %>% select(Tx_name, n_SNMs, text_pos)
text_backing_2$x_min <- text_backing_2$text_pos - 2.5
text_backing_2$x_max <- text_backing_2$text_pos + 2.5
text_backing_2$y_min <- text_backing_2$n_SNMs - 0.45
text_backing_2$y_max <- text_backing_2$n_SNMs + 0.45

SNM_n_rect_plot <- SNMs_final_counts_table %>%
  ggplot() + 
  # geom_hline(data = SNMs_final_counts_stats, aes(yintercept = mean_n, color = Tx_name), size = 0.75) +
  geom_rect(data = text_backing, aes(xmin = x_min, xmax = x_max, ymin = y_min, ymax = y_max), fill = "white") +
  geom_rect(aes(xmin = x_left, xmax = x_right, ymin = y_bot, ymax = y_top, fill = Tx_name)) +
  geom_crossbar(data = SNMs_means, aes(x = 0, y = m, ymin = m, ymax = m), size = 0.1, width = n_muts_lim) +
  geom_rect(data = text_backing_2, aes(xmin = x_min, xmax = x_max, ymin = y_min, ymax = y_max), fill = txPal[2]) +
  geom_text(aes(x = text_pos, y = n_SNMs, 
                 label = n_muts), color = SNMs_final_counts_table$text_color, size = 6) +
  scale_x_continuous(breaks = 0, label = NULL, limits = c(-n_muts_lim/2 - 2, n_muts_lim/2 + 2)) +
  scale_y_continuous(breaks = seq.int(0, max(SNMs_final_counts_table$n_SNMs, 1)), 
                     name = "Number of Point Mutations") +
  scale_fill_manual(values = txPal) + 
  scale_color_manual(values = txPal) +
  # coord_flip() +
  theme(legend.position = "none",
        panel.grid.minor.y = element_blank(),
        plot.margin = unit(c(t = 5, r = 5, b = 5, l = 8), "mm"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(vjust = 3),
        strip.text = element_blank(),
        # strip.text.x = element_blank(),
        text = element_text(size = 20)) +
  facet_grid(.~Tx_name)

SNM_n_rect_plot

ggsave(file.path(outIntDir, "SNMcount_freq_2022_02.png"), 
       plot = SNM_n_rect_plot,
       device = "png",
       width = 16, height = 9, 
       units = "in",
       dpi = 600)

# Segmented center-out vertical histogram
###############################################################################
SNMs_final_counts_noZero <- SNMs_final_counts_table %>% filter(n_muts > 0)
n_ea_SNM <- SNMs_final_counts_noZero %>% pull(n_muts)
clone_i <- unlist(sapply(n_clones_xTx$n, function(x) 1:x))
clone_Tx <- unlist(mapply(rep, times = n_clones_xTx$n, x = Tx_name_levels))
clone_x <- unlist(sapply(n_ea_SNM, function(x) 1:x))
clone_y <- unlist(mapply(rep, times = SNMs_final_counts_table$n_muts, 
                         x = SNMs_final_counts_table$n))
clone_x_adj <- SNMs_final_counts_noZero %>% 
  group_by(Tx_name, n) %>% 
  summarize(clone_i = 1:n_muts, adj = 1:n_muts - (max(n_muts) + 1) / 2)

SNMs_counts_plot <- data.frame(Tx_name = clone_Tx, clone = clone_i, 
                               clone_xi = clone_x, clone_n_SNMs = clone_y,
                               x_adj = clone_x_adj$adj)
SNMs_counts_plot$Tx_name <- factor(SNMs_counts_plot$Tx_name, levels = Tx_name_levels)


SNM_n_tile_plot <- SNMs_counts_plot %>%
  ggplot() + 
  # geom_hline(data = SNMs_final_counts_stats, aes(yintercept = mean_n, color = Tx_name), size = 0.75) +
  geom_tile(aes(x = x_adj, y = clone_n_SNMs, fill = Tx_name), width = 0.8, height = 0.8) +
  # geom_label(data = SNMs_final_counts_table %>% filter(n_muts > 0),
  #            aes(x = x_center, y = n_SNMs, label = n_muts), 
  #            size = 6, label.size = 0) +
  scale_x_continuous(breaks = 0, label = "") +
  scale_y_continuous(breaks = seq.int(0, max(SNMs_final_counts_table$n, 1))) +
  # xlab("Number of Clones") + 
  ylab("Number of Point Mutations") +
  scale_fill_manual(values = txPal, name = "Treatment") + 
  scale_color_manual(values = txPal, name = "Treatment") +
  # coord_flip() +
  theme(legend.position = "none",
        panel.grid.minor.y = element_blank(),
        axis.title.x = element_blank(),
        # strip.text.x = element_blank(),
        text = element_text(size = 14)) +
  facet_grid(.~Tx_name, switch = "x")

SNM_n_tile_plot

ggsave(file.path(outIntDir, "SNMcount_freq_2022_02.png"), 
       plot = SNM_n_plot,
       device = "png",
       width = 11, height = 8.5, 
       units = "in",
       dpi = 600)


SNM_n_plot_line <- SNMs_final_counts_table %>%
  ggplot() + 
  geom_vline(data = SNMs_means, aes(xintercept = m, color = Tx_name), show.legend = F) +
  geom_label_repel(data = SNMs_means, aes(x = m, y = 0, label = round(m, 2), color = Tx_name), 
                   hjust = "left", label.size = 0, size = 3.5, show.legend = F) +
  geom_line(aes(x = n_SNMs, y = density, color = Tx_name)) +
  geom_point(aes(x = n_SNMs, y = density, color = Tx_name)) +
  # geom_label(data = lab_df, aes(x = n_SNMs, y = n_muts), label = SNM_counts_label,
  #            hjust = "left", label.size = 0, size = 4) +
  # scale_x_continuous(breaks = seq(0, 15, 1)) +
  xlab("Number of Single Nucleotide Mutations") + ylab("Proportion of Clones") +
  scale_fill_manual(values = txPal, name = "Treatment") + 
  scale_color_manual(values = txPal, name = "Treatment") +
  theme(legend.position = "bottom", text = element_text(size = 20))

SNM_n_plot_line

ggsave(file.path(outIntDir, "SNMcount_freq_2021_12_v3.png"), 
       plot = SNM_n_plot,
       device = "png",
       width = 11, height = 8.5, 
       units = "in",
       dpi = 600)

SNM_rate_plot <- SNMs_final_counts_stats %>%
  ggplot(aes(x = Tx_name, y = rate)) + 
  # geom_vline(data = SNMs_means, aes(xintercept = m, color = Tx_name), show.legend = F) +
  # geom_label_repel(data = SNMs_means, aes(x = m, y = 0, label = round(m, 2), color = Tx_name), 
  #                  hjust = "left", label.size = 0, size = 3.5, show.legend = F) +
  geom_errorbar(aes(ymin = rateCI95lo, ymax = rateCI95hi, color = Tx_name)) +
  geom_point(aes(color = Tx_name)) +
  # geom_label(data = lab_df, aes(x = n_SNMs, y = n_muts), label = SNM_counts_label,
  #            hjust = "left", label.size = 0, size = 4) +
  # scale_x_continuous(breaks = seq(0, 15, 1)) +
  ylim(0,NA) +
  xlab("Strain") + ylab("Point mutation rate (+/- 95%CI)") +
  scale_fill_manual(values = txPal, name = "Treatment") + 
  scale_color_manual(values = txPal, name = "Treatment") +
  theme(legend.position = "bottom", text = element_text(size = 20))

SNM_rate_plot

ggsave(file.path(outIntDir, "SNMrate_2021_12_v3.png"), 
       plot = SNM_rate_plot,
       device = "png",
       width = 11, height = 8.5, 
       units = "in",
       dpi = 600)



###############################################################################
# Do the SNM distributions follow a Poisson?

# Poisson test for pooled data
nSNMs_freq <- SNMs_final_counts %>% count(n_SNMs = n)
mean_nSNMs <- mean(SNMs_final_counts$n)
pois_dist <- dpois(nSNMs_freq$n_SNMs, mean_nSNMs, log = F)
resid <- 1-sum(pois_dist)
resid_prop <- (pois_dist+resid/length(pois_dist))*resid
pois_dist_crct <- pois_dist + resid_prop

chi_pois_test <- chisq.test(nSNMs_freq$n, p=pois_dist_crct)
chi_posi_label <- paste0("Chi Sq. p = ", ifelse(chi_pois_test$p.value < 1E-15, "2.2e-16", 
                                                formatC(chi_pois_test$p.value, 2, format = "E")))


ggplot() + geom_col(aes(x=nSNMs_freq$n_SNMs, y=nSNMs_freq$n/sum(nSNMs_freq$n))) +
  geom_line(aes(x=nSNMs_freq$n_SNMs, y=pois_dist_crct), color="red") +
  geom_label(aes(x = 4, y = 0.3, label = chi_posi_label), label.size = 0) +
  xlab("Number of SNMs") + ylab("Fraction of clones")

#  Poisson test for Tx data
# Generate Poisson dist for plotting 
pois_dist_tx_crct <- list()
for (t in 1:3) {
  # t = 1
  tx <- c("WT", "Cas9", "Drive")[t]
  idx <- SNMs_final_counts_table$Tx_name == tx
  idx2 <- SNMs_final_counts$Tx_name == tx
  pois_dist_tx <- dpois(0:max(SNMs_final_counts_table$n_muts), mean(SNMs_final_counts$n[idx2]), log = F)
  resid <- 1-sum(pois_dist_tx)
  resid_prop <- (pois_dist_tx+resid/length(pois_dist_tx))*resid
  pois_dist_tx_crct[[t]] <- pois_dist_tx + resid_prop
}

data_len <- length(pois_dist_tx_crct[[1]])

pois_df <- data.frame(Tx_name = c(rep("WT", data_len), 
                                  rep("Cas9", data_len), 
                                  rep("Drive",data_len)), 
                      value = rep(0:(data_len - 1), 3), 
                      prop = unlist(pois_dist_tx_crct))

pois_df$Tx_name <- factor(pois_df$Tx_name, levels = c("WT", "Cas9", "Drive"))

# Estimate p-values for null hyp of equal variances by Tx ------
nTrials <- 10000
var_list <- list()
var_sample <- list()
p_vals <- data.frame(NULL)
for (t in 1:3) {
  tx <- c("WT", "Cas9", "Drive")[t]
  idx2 <- SNMs_final_counts$Tx_name == tx
  var_sample[t] <- var(SNMs_final_counts$n[idx2])
  var_pois <- c()
  for (r in 1:nTrials) {
    pois_sample <- rpois(length(SNMs_final_counts$n[idx2]), 
                         mean(SNMs_final_counts$n[idx2]))
    var_pois[r] <- var(pois_sample)
  }
  var_list[[t]] <- var_pois
  p_vals <- rbind(p_vals, c(Tx_name = tx, p = signif(sum(var_list[[t]] > var_sample[[t]])/nTrials, 3)))
  # p_vals <- c(p_vals, sum(var_list[[t]] > var_sample[[t]])/nTrials)
  # names(p_vals[[t]]) <- tx
}

colnames(p_vals) <- c("Tx_name", "p")
p_vals$Tx_name <- factor(p_vals$Tx_name, levels = c("WT", "Cas9", "Drive"))
p_vals$varLabels <- c(paste0(ifelse(p_vals$p != 0,  "var p = ",  "var p < "), ifelse(p_vals$p != 0, p_vals$p, 1/nTrials)))

Chi_pois_tx <- list()
for (t in 1:3) {
  # t=3
  tx <- levels(SNMs_final_counts_table$Tx_name)[t]
  idx <- SNMs_final_counts_table$Tx_name == tx
  idx2 <- SNMs_final_counts$Tx_name == tx
  pois_dist <- dpois(SNMs_final_counts_table$n_muts[idx], mean(SNMs_final_counts$n[idx2]), log = F)
  resid <- 1-sum(pois_dist)
  resid_prop <- (pois_dist+resid/length(pois_dist))*resid
  pois_dist_crct <- pois_dist + resid_prop
  Chi_pois_tx[[t]] <- chisq.test(x = SNMs_final_counts_table$n_muts[idx], p = pois_dist_crct)
}

# Chi_pois_tx[[1]]$p.value
p_vals$chiLabels <- c(paste0("X-sq p = ", 
                             signif(c(Chi_pois_tx[[1]]$p.value, 
                                      Chi_pois_tx[[2]]$p.value, 
                                      Chi_pois_tx[[3]]$p.value), 3)))
# p_vals$pLabels <- paste0()

nSNMs_isPois <- SNMs_final_counts_table %>% 
  ggplot() + geom_col(aes(x = n, y = density)) +
  geom_line(data = pois_df, aes(x = value, y = prop), color="red", inherit.aes = F) +
  geom_label(data = p_vals, x = data_len * 0.75, y = 0.35, label.size = 0,
             size = 4, label.padding = unit(0.5, "lines"),
             aes(label = paste0(chiLabels, "\n", varLabels))) +
  theme(strip.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12)) +
  #annotate(geom = "text", x = 25, y = 0.10, label=p_vals$p) +
  xlab("# LOH events") + ylab("Proportion") +
  facet_grid(Tx_name~.)

nSNMs_isPois


final_mut_set %>% ggplot() + 
  # scale_x_continuous(labels=as.character(chrom_bound_BY$CHROM),
  #                    breaks=(chrom_bound_BY$Start+chrom_bound_BY$End)/2) +
  # annotate(geom="rect", xmin=chrom_bound_BY$Start[c(TRUE, FALSE)],
  #          xmax=chrom_bound_BY$End[c(TRUE, FALSE)], ymin=-Inf, ymax=Inf, alpha=0.05, color = "grey90") +
  geom_jitter(aes(x = POS, y = Tx_name), height = 0.1, width = 0) + 
  facet_wrap(~CHROM, ncol = 4, scales = "free_x") +
  theme(panel.background = element_rect(color = "grey50"))



# Point mutations along genome !!!!!!!!########################################
chrom_lengths_df <- data.frame(CHROM = chrom_bound_BY$CHROM,
                               Start = 0, End = chrom_lengths_BY/1000)
final_mut_set$POS_kb <- final_mut_set$POS/1000
chrom_lim <- max(chrom_lengths_BY)/1000
SNMs_pos_plot <- final_mut_set %>% ggplot() +
  # scale_x_continuous(labels=as.character(chrom_bound_BY$CHROM),
  #                    breaks=(chrom_bound_BY$Start+chrom_bound_BY$End)/2,
  #                    expand = c(0,0)) +
  # annotate(geom="rect", xmin=chrom_lengths_df$Start,
  #          xmax=chrom_lengths_df$End, ymin=-Inf, ymax=Inf, alpha=0.05, color = "grey90") +
  geom_segment(data = chrom_lengths_df, aes(x = Start, xend = End, y = CHROM, yend = CHROM), color = "grey20", alpha = 0.5) +
  geom_point(aes(x = POS_kb, y = CHROM, color = Tx_name, fill = Tx_name), shape = 21, size = 4, alpha = 0.2) + 
  geom_point(aes(x = POS_kb, y = CHROM, color = Tx_name, fill = Tx_name), shape = 1, size = 4) + 
  scale_x_continuous(limits = c(0, chrom_lim), 
                     breaks = seq(0, chrom_lim, 100), 
                     expand = c(0, 50), name = "Position (kbp)") +
  scale_color_manual(values = txPal) +
  scale_fill_manual(values = txPal) +
  scale_y_discrete(limits = rev(levels(final_mut_set$CHROM)), name = "Chromosome") +
  # facet_wrap(~CHROM, ncol = 1) +
  theme(legend.position = "bottom",
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    # panel.background = element_rect(color="grey80"),
    text = element_text(size = 14))

SNMs_pos_plot

ggsave(file.path(outIntDir, "SNMs_pos_plot_2022_02.png"), 
       plot = SNMs_pos_plot,
       device = "png",
       width = 11, height = 8.5, 
       units = "in",
       dpi = 600)



# UNUSED ######################################################################

mut_ID_CHROM_count <- final_mut_set %>% 
  group_by(ID) %>% 
  count(CHROM, name = "n_SNMs", .drop = F) %>% arrange(ID, CHROM) %>% as.data.frame()

mut_ID_CHROM_filter <- mut_ID_CHROM_count %>% summarize(ID_CHROM = paste0(ID, "_", CHROM))
mut_ID_CHROM_exclude <- mut_ID_CHROM_filter$ID_CHROM %in% low_cover_ID_CHROM$ID_CHROM
mut_ID_CHROM_count <- mut_ID_CHROM_count %>% filter(!mut_ID_CHROM_exclude)
mut_ID_CHROM_count <- CategoriesFromID(mut_ID_CHROM_count)

Tx_mut_CHROM_count <- mut_ID_CHROM_count %>%
  # dplyr::ungroup() +
  group_by(Tx_name, CHROM) %>%
  dplyr::summarize(clone_rate = mean(n_SNMs), se_rate = se(n_SNMs),
                   tot_SNMs = sum(n_SNMs), .groups = "keep") 

Tx_mut_CHROM_count$gen_rate <- Tx_mut_CHROM_count$clone_rate / n_gens
Tx_mut_CHROM_count$bp_rate <- Tx_mut_CHROM_count$gen_rate / rep(chrom_lengths_BY, times = 3)
Tx_mut_CHROM_count$gen_se <- Tx_mut_CHROM_count$se_rate / n_gens
Tx_mut_CHROM_count$bp_se <- Tx_mut_CHROM_count$gen_se / rep(chrom_lengths_BY, times = 3)
Tx_mut_CHROM_count$bp_up <- Tx_mut_CHROM_count$bp_rate + Tx_mut_CHROM_count$bp_se * 1.96
Tx_mut_CHROM_count$bp_down <- ifelse(Tx_mut_CHROM_count$bp_rate > Tx_mut_CHROM_count$bp_se * 1.96,
                                     Tx_mut_CHROM_count$bp_rate - Tx_mut_CHROM_count$bp_se * 1.96,
                                     0)

# Test difference in chromosomal distributions of mutations between treatment pairs

WT_chrom_count <- Tx_mut_CHROM_count %>% filter(Tx_name == "WT") %>% pull(tot_SNMs)
Cas9_chrom_count <- Tx_mut_CHROM_count %>% filter(Tx_name == "Cas9") %>% pull(tot_SNMs)
Drive_chrom_count <- Tx_mut_CHROM_count %>% filter(Tx_name == "Drive") %>% pull(tot_SNMs)
chisq.test(WT_chrom_count, Cas9_chrom_count)
chisq.test(WT_chrom_count, Drive_chrom_count)
chisq.test(Cas9_chrom_count, Drive_chrom_count)

# Test significance of treatment pairs by chromosome
Tx_combo <- data.frame(Tx_1 = c("WT", "WT", "Cas9"), Tx_2 = c("Cas9", "Drive", "Drive"))

mut_CHROMrate_perm <- data.frame(NULL)
for(tx in 1:nrow(Tx_combo)){
  # tx = 1
  all_Chr_perm <- data.frame(Combo = with(Tx_combo[tx, ], paste(Tx_1, Tx_2, sep = "-")), 
                             CHROM = chrom_bound_BY$CHROM, 
                             obsvStat = 0, critVal = 0, pVal = 0, rejectNull = 0)
  
  for(c in seq_along(chrom_bound_BY$CHROM)) {
    # c = 3
    chr = chrom_bound_BY$CHROM[c]
    Chr_SNM_counts_ID <- mut_ID_CHROM_count %>% filter(CHROM == chr)
    Chr_perm <- perm_test(Chr_SNM_counts_ID, 
                          cat_var = "Tx_name", 
                          cat_names = with(Tx_combo[tx, ], c(Tx_1, Tx_2)), 
                          response_var = "n_SNMs", 
                          n_perms = 10000, alt_hyp = "two-tailed")
    all_Chr_perm[c, 3:6] <- unlist(Chr_perm[1:4])
  }
  mut_CHROMrate_perm <- rbind(mut_CHROMrate_perm, all_Chr_perm)
}

mut_CHROMrate_perm <- BHcorrection(mut_CHROMrate_perm)

mut_CHROMrate_perm$sig_lab_p <- ifelse(mut_CHROMrate_perm$rejectNull == 1, "*", "")
mut_CHROMrate_perm$sig_lab_BH <- ifelse(mut_CHROMrate_perm$rejectNull_BH == 1, "**", "")
mut_CHROMrate_perm$sig_lab_both <- paste0(mut_CHROMrate_perm$sig_lab_p, 
                                          ifelse(mut_CHROMrate_perm$rejectNull_BH == 1, ", ", ""),
                                          mut_CHROMrate_perm$sig_lab_BH)

combo_cols <- colsplit(mut_CHROMrate_perm$Combo, "-", c("Tx1", "Tx2"))
mut_CHROMrate_perm <- cbind(mut_CHROMrate_perm, combo_cols)
mut_CHROMrate_perm$Tx1 <- factor(mut_CHROMrate_perm$Tx1, levels = c("WT", "Cas9", "Drive"))
mut_CHROMrate_perm$Tx2 <- factor(mut_CHROMrate_perm$Tx2, levels = c("WT", "Cas9", "Drive"))
mut_CHROMrate_perm$Combo <- factor(mut_CHROMrate_perm$Combo)
mut_sig_CHROMrate <- mut_CHROMrate_perm %>% filter(rejectNull == 1) 
Tx_mut_CHROM_count$rom_CHROM <- chrom_IDs$rom_CHROM[as.numeric(Tx_mut_CHROM_count$CHROM)]
y_lab <- max(Tx_mut_CHROM_count$bp_up) * 1.05

DNrate_CHROM_line <- Tx_mut_CHROM_count %>% 
  ggplot() + 
  geom_errorbar(aes(x = rom_CHROM + , ymin = bp_rate, ymax = bp_rate, color = Tx_name),
                width = 0.6, position = position_dodge(width = 0.6), size = 0.75) +
  geom_point(aes(x = ))
  geom_errorbar(aes(x = Tx_name, ymin = bp_down, ymax = bp_up, color = Tx_name),
                width = 0, position = position_dodge(width = 0.6), size = 0.75) +
  # geom_bracket(data = mut_sig_CHROMrate, 
  #              y.position = y_lab + (as.numeric(mut_sig_CHROMrate$Combo)-1)*1E-11, step.increase = 0, 
  #              vjust = 0.4, label.size = 5,
  #              aes(xmin = Tx1, xmax = Tx2, label = sig_lab_both)) +
  geom_text(aes(x = Tx_name, y = -1E-11, label = tot_SNMs, color = Tx_name), 
            position = position_dodge(width = 0.6), show.legend = F) +
  scale_color_manual(values = txPal, name = "Drive Type") +
  scale_fill_manual(values = txPal, name = "Drive Type", guide = "none") +
  # ylim(c(0, NA)) +
  ylab("SNM mutation rate (/bp/generation)") + 
  xlab("Chromosome") +
  facet_wrap(~CHROM, ncol = 16, strip.position = "bottom") +
  # xlab("Chromosome")  +
  theme(axis.text.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.background = element_rect(color = "grey50"),
        text = element_text(size = 14))

DNrate_CHROM_line

ggsave(file.path(outIntDir, "SNMrate_CHROM_line_2021_11.png"), 
       plot = DNrate_CHROM_line,
       device = "png",
       width = 16, height = 9, 
       units = "in",
       dpi = 600)

# Compare SNMs by CHROM across Tx_name
final_mut_CHROM_count <- final_mut_set %>% 
  ungroup() %>% 
  count(Tx_name, CHROM) %>% 
  arrange(Tx_name, CHROM)

final_mut_CHROM_count_wide <- final_mut_CHROM_count %>% 
  pivot_wider(id_cols = CHROM, names_from = Tx_name, values_from = n) %>%
  as.data.frame()
final_mut_CHROM_count_wide[is.na(final_mut_CHROM_count_wide)] <- 0
final_mut_CHROM_count_wide <- final_mut_CHROM_count_wide %>%
  arrange(CHROM)

final_mut_CHROM_count <- final_mut_CHROM_count_wide %>% 
  pivot_longer(cols = -CHROM, names_to = "Tx_name", values_to = "n_SNMs") 
final_mut_CHROM_count$Tx_name <- factor(final_mut_CHROM_count$Tx_name, 
                                        levels = c("WT", "Cas9", "Drive"))
final_mut_CHROM_count <- final_mut_CHROM_count %>%
  arrange(Tx_name, CHROM)

final_mut_CHROM_count <- final_mut_CHROM_count %>% 
  group_by(Tx_name) %>%
  mutate(density = n_SNMs/sum(n_SNMs)) %>% ungroup() %>%
  as.data.frame()

final_mut_CHROM_count$clone_rate <- final_mut_CHROM_count$n_SNMs / rep(n_clones_xTx$n, each = 16)
final_mut_CHROM_count$gen_rate <- final_mut_CHROM_count$clone_rate / n_gens

final_mut_CHROM_count$bp_rate <- final_mut_CHROM_count$gen_rate / rep(chrom_lengths_BY, times = 3)

aov_chom_rate <- aov(bp_rate ~ Tx_name + CHROM, data = final_mut_CHROM_count)
aov_Tx_p <- unlist(summary(aov_chom_rate))["Pr(>F)1"]
aov_chrom_p <- unlist(summary(aov_chom_rate))["Pr(>F)2"]

aov_chrom_label <- paste0("ANOVA \n  Drive Type p = ", signif(aov_Tx_p, 3),
                          "\n  Chromosome p = ", signif(aov_chrom_p, 3))

N_chrom_count <- final_mut_CHROM_count_wide[, "WT"]
H_chrom_count <- final_mut_CHROM_count_wide[, "Cas9"]
F_chrom_count <- final_mut_CHROM_count_wide[, "Drive"]

chi_chrom_N_H <- chisq.test(cbind(N_chrom_count, H_chrom_count))
chi_chrom_N_F <- chisq.test(cbind(N_chrom_count, F_chrom_count))
chi_chrom_H_F <- chisq.test(cbind(H_chrom_count, F_chrom_count))


chi_chrom_label <- paste0("Chi Sq. \n  WTvCas9 p = ", signif(chi_chrom_N_H$p.value, 3),
                          "\n  WTvDrive p = ", signif(chi_chrom_N_F$p.value, 3),
                          "\n  Cas9vDrive p = ", signif(chi_chrom_H_F$p.value, 3))

both_label <- paste0(aov_chrom_label, "\n", chi_chrom_label)

lab_y <- 0.9*max(final_mut_CHROM_count$gen_rate)
final_mut_CHROM_count %>% 
  ggplot() + 
  geom_line(aes(x = CHROM, y = gen_rate, color = Tx_name, group = Tx_name),
            size = 0.15, alpha = 0.5) +
  geom_point(aes(x = CHROM, y = gen_rate, color = Tx_name), 
             size = 3) +
  # geom_label(aes(x = 14, y = lab_y), label = both_label, hjust = 0) +
  scale_color_manual(values = txPal) +
  xlab("Chromosome") + ylab("SNM rate /chrom/clone") +
  theme(panel.grid.minor.y = element_blank())


###############################################################################
# Compare SNMs across genome
ecdf_in <- final_mut_set

ad_test_WTvCas9 <- kSamples::ad.test(ecdf_in %>% filter(Tx_name == "WT") %>% pull(POSi), 
                                     ecdf_in %>% filter(Tx_name == "Cas9") %>% pull(POSi))

ad_test_WTvDrive <- kSamples::ad.test(ecdf_in %>% filter(Tx_name == "WT") %>% pull(POSi), 
                                      ecdf_in %>% filter(Tx_name == "Drive") %>% pull(POSi))

ad_test_Cas9vDrive <- kSamples::ad.test(ecdf_in %>% filter(Tx_name == "Cas9") %>% pull(POSi), 
                                        ecdf_in %>% filter(Tx_name == "Drive") %>% pull(POSi))

ad_label <- paste0("AD test \n  WTvCas9 p = ", round(max(ad_test_WTvCas9$ad[, 3]), 3),
                   "\n  WTvDrive p = ", round(max(ad_test_WTvDrive$ad[, 3]), 3),
                   "\n  Cas9vDrive p = ", round(max(ad_test_Cas9vDrive$ad[, 3]), 3))

ecdf_in %>% ggplot() +
  # scale_x_continuous(labels=as.character(chrom_bound_BY$CHROM),
  #                    breaks=(chrom_bound_BY$Start+chrom_bound_BY$End)/2,
  #                    expand = c(0,0)) +
  # annotate(geom="rect", xmin=chrom_bound_BY$Start[c(TRUE, FALSE)],
  #          xmax=chrom_bound_BY$End[c(TRUE, FALSE)], ymin=-Inf, ymax=Inf, alpha=0.05, color = "grey90") +
  # annotate(geom="segment", x = 0, y = 0, xend = chrom_bound_BY$End[16], yend = 1, size = 0.25) +
  # geom_segment(aes(x = 0, y = 0, xend = chrom_bound_BY$End[16], yend = 1), size = 0.1, alpha = 0.5) +
  # geom_label(aes(x = chrom_bound_BY$End[1], y = 0.85, label = ad_label), hjust = "left") +
  facet_wrap(~CHROM, ncol = 4, scales = "free_x") +
  stat_ecdf(aes(x = POS, color = Tx_name)) + 
  # geom_text(data = Tx_mut_CHROM_count, 
  #           aes(x = 10000, y = 1.05 - as.numeric(Tx_name)*0.1, label = tot_SNMs, color = Tx_name), 
  #           position = position_dodge(width = 0.6), show.legend = F) +
  # stat_ecdf(aes(x = POSi), color = "red", size = 0.3) + 
  scale_color_manual(values = txPal) + 
  xlab("Genome Position") + ylab("Cumulative Fraction of SNM mutations") +
  theme(panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank(),
        panel.background = element_rect(color="grey80"),
        text = element_text(size = 14))


library(goftest)
ad_test_unif <- goftest::ad.test(ecdf_in %>% filter(Tx_name == "WT") %>% pull(POSi), 
                                 null = "punif")
ad_label_unif <- paste0("AD test \n  Uniform dist. p = ", formatC(ad_test_unif$p.value, format = "E", digits = 2))

ecdf_in %>% ggplot() +
  scale_x_continuous(labels=as.character(chrom_bound_BY$CHROM),
                     breaks=(chrom_bound_BY$Start+chrom_bound_BY$End)/2,
                     expand = c(0,0)) +
  annotate(geom="rect", xmin=chrom_bound_BY$Start[c(TRUE, FALSE)],
           xmax=chrom_bound_BY$End[c(TRUE, FALSE)], ymin=-Inf, ymax=Inf, alpha=0.05, color = "grey90") +
  geom_segment(aes(x = 0, y = 0, xend = chrom_bound_BY$End[16], yend = 1), size = 0.1) +
  geom_label(aes(x = chrom_bound_BY$End[1], y = 0.85, label = ad_label_unif), hjust = "left") +
  stat_ecdf(aes(x = POSi), color = "red", size = 0.3) + 
  scale_color_manual(values = txPal) + 
  xlab("Genome Position") + ylab("Cumulative Fraction of SNM mutations") +
  theme(panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank(),
        panel.background = element_rect(color="grey80"),
        text = element_text(size = 14))

chrom_lengths_BY_df <- data.frame(CHROM = unique(denovo_SNMs$CHROM), length = chrom_lengths_BY)

ecdf_in %>% 
  # arrange(Tx_name, CHROM, POS) %>%
  # group_by(Tx_name, CHROM) %>% 
  arrange(Tx_name, POSi) %>%
  group_by(Tx_name) %>% 
  mutate(cum_n = row_number()) %>%
  ggplot() +
  # facet_wrap(~CHROM, ncol = 4, scales = "free_x") +
  # geom_segment(data = chrom_lengths_BY_df, aes(x = 0, y = 0, xend = length, yend = 1), size = 0.25) +
  # stat_ecdf(aes(x = POS), color = "red", size = 0.3) +
  geom_step(aes(x = POSi, y = cum_n, color = Tx_name)) +
  # stat_ecdf(aes(x = POS, color = Tx_name), size = 0.3) +
  scale_x_continuous(expand = c(0,0)) +
  scale_color_manual(values = txPal) + 
  xlab("Chromosome Position") + 
  ylab("Cumulative Fraction of SNM mutations") +
  theme(
    # panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.background = element_rect(color="grey80"),
    text = element_text(size = 14))


