# LOH position and length distribution analysis


# Distribution of LOH lengths by treatment and interstitial or terminal
###############################################################################
all_LOHbounds_merge_NS$Tx_ID <- Recode_Tx_ID(all_LOHbounds_merge_NS$Tx, "Tx")
Length_data_in <- all_LOHbounds_merge_NS

Length_data_in$Tx_ID <- Recode_Tx_ID(Length_data_in$Tx)

Length_data_in$type <- factor(ifelse(Length_data_in$isTerm, 
                                          "Terminal", "Interstitial"))

LOH_type_xTx <- all_LOHcounts_merge_NS %>% 
  group_by(Tx_ID) %>% 
  summarize(n_iLOH = sum(n_Inter),
            n_tLOH = sum(n_Term),
            f_iLOH = sum(n_Inter)/sum(n_LOH),
            mean_iLOH = mean(n_Inter),
            mean_tLOH = mean(n_Term))

LOH_type_xID <- all_LOHcounts_merge_NS %>% 
  group_by(Tx_ID, ID) %>% 
  summarize(n_iLOH = sum(n_Inter),
            n_tLOH = sum(n_Term),
            f_iLOH = sum(n_Inter)/sum(n_LOH))

  
iLOH_CIs <- data.frame(NULL)
for(tx in Tx_ID_levels) {
  tx_n_LOH <- LOH_type_xID %>% filter(Tx_ID == tx) %>% pull(n_iLOH)
  tx_boot <- boot(tx_n_LOH, statistic = mean.fun, R = 10000)
  tx_boot_CI <- boot.ci(boot.out = tx_boot, conf = 0.95,
                        type = "basic")$basic[4:5]
  tx_CI <- data.frame(Tx_ID = tx, low95CI = tx_boot_CI[1], up95CI = tx_boot_CI[2])
  iLOH_CIs <- rbind(iLOH_CIs, tx_CI)
}

tLOH_CIs <- data.frame(NULL)
for(tx in Tx_ID_levels) {
  tx_n_LOH <- LOH_type_xID %>% filter(Tx_ID == tx) %>% pull(n_tLOH)
  tx_boot <- boot(tx_n_LOH, statistic = mean.fun, R = 10000)
  tx_boot_CI <- boot.ci(boot.out = tx_boot, conf = 0.95,
                        type = "basic")$basic[4:5]
  tx_CI <- data.frame(Tx_ID = tx, low95CI = tx_boot_CI[1], up95CI = tx_boot_CI[2])
  tLOH_CIs <- rbind(tLOH_CIs, tx_CI)
}

LOH_type_xTx <- merge(LOH_type_xTx, iLOH_CIs, by = "Tx_ID")
LOH_type_xTx <- merge(LOH_type_xTx, tLOH_CIs, by = "Tx_ID", suffixes = c("_i", "_t"))
LOH_type_xTx <- LOH_type_xTx %>% arrange(Tx_ID)
LOH_type_xID$Tx_ID <- Recode_Tx_ID(LOH_type_xID$Tx_ID, "Tx_ID")

iLOH_xn_count <- LOH_type_xID %>% count(Tx_ID, n_iLOH)
tLOH_xn_count <- LOH_type_xID %>% count(Tx_ID, n_tLOH)

iLOH_xn_count$d_clones <- operate_by_factor_match(n_clones_LOH_xTx[, c("Tx_ID", "n")], 
                                                iLOH_xn_count[, c("Tx_ID", "n")], 
                                               .fun = function(x, y) y/x)

tLOH_xn_count$d_clones <- operate_by_factor_match(n_clones_LOH_xTx[, c("Tx_ID", "n")], 
                                                  tLOH_xn_count[, c("Tx_ID", "n")], 
                                                  .fun = function(x, y) y/x)

chisq.test(t(LOH_type_xTx[, 2:3]))

LOH_type_aov <- aov(f_iLOH ~ Tx_ID, data = LOH_type_xID)
summary(LOH_type_aov)

Tx_combos <- expand.grid(Tx_1 = c("W", "C", "D"), 
                         Tx_2 = c("W", "C", "D"), 
                         LOH = c("n_iLOH", "n_tLOH")) %>% 
  filter(as.numeric(Tx_1) < as.numeric(Tx_2)) %>% 
  mutate(Tx_1 = as.character(Tx_1),
         Tx_2 = as.character(Tx_2),
         LOH = as.character(LOH),
         combo = paste0(Tx_1, "_", Tx_2, "_", LOH))

Tx_combos <- expand.grid(Tx_1 = c("W"), 
                         Tx_2 = c("C", "D"), 
                         LOH = c("n_iLOH", "n_tLOH")) %>% 
  # filter(as.numeric(Tx_1) < as.numeric(Tx_2)) %>% 
  mutate(Tx_1 = as.character(Tx_1),
         Tx_2 = as.character(Tx_2),
         LOH = as.character(LOH),
         combo = paste0(Tx_1, "_", Tx_2, "_", LOH))

Tx_combos_list <- nest(.data = Tx_combos, combos = -combo)$combos %>% as.list()
names(Tx_combos_list) <- Tx_combos$combo

LOH_type_p <- lapply(Tx_combos_list, function(x) perm_test(LOH_type_xID,
                                             cat_var = "Tx_ID", 
                                             cat_names = c(unlist(x[1]), unlist(x[2])),
                                             response_var = unlist(x[3]),
                                             rtrn = "p"))

p.adjust(unlist(LOH_type_p), method = "BH")
# LOH_type_p_df <- data.frame(combo = names(Tx_combos_list), 
#                             p_value = matrix(LOH_type_p))
# BHcorrection(LOH_type_p_df, "p_value")
LOH_type_xID %>% 
  pivot_longer(cols = c(n_iLOH, n_tLOH), 
               names_to = "type", values_to = "n_LOH") %>%
  mutate(type = factor(type)) %>%
  group_by(Tx_ID) %>%
  group_modify(~ as.data.frame(lm(n_LOH ~ type, data = .x)))

summary(lm(n_tLOH ~ n_iLOH + Tx_ID, data = LOH_type_xID %>% filter(Tx_ID != "C")))

LOH_type_plot <- LOH_type_xID %>% 
  filter(Tx_ID != "D") %>%
  ggplot(aes(x = n_iLOH, y = n_tLOH, color = Tx_ID)) + 
  stat_smooth(method = "lm", fill = "grey70", alpha = 0.2) +
  geom_jitter(height = 0.3, width = 0.3, size = 3.5, alpha = 0.8) +
  scale_x_continuous(breaks = 0:10, name = "Number of interstitial LOH events") +
  scale_y_continuous(breaks = 0:10, name = "Number of terminal LOH events") +
  # xlab("Number of interstitial LOH events") +
  # ylab("Number of terminal LOH events") +
  scale_color_manual(values = txPal, name = "Strain") +
  theme(panel.grid.minor = element_blank(),
        text = element_text(size = 28),
        legend.position = c(0.1, 0.9),
        legend.background = element_rect(fill = "white", color = "white"))

LOH_type_plot

y_pos_mean_i <- max(iLOH_xn_count$n)
y_space_i <- 2

iLOH_density <- iLOH_xn_count %>% 
  filter(Tx_ID != "D") %>%
  ggplot() +
  geom_col(aes(x = n_iLOH, y = n, fill = Tx_ID, color = Tx_ID), 
                 alpha = 0.7,
                 position = position_dodge(preserve = "single", width = 0.8), width = 0.8) +
  # geom_point(data = LOH_type_xTx, 
  #            aes(x = mean_iLOH, y = y_pos_mean_i + (4 - as.numeric(Tx_ID)) * y_space_i,
  #                color = Tx_ID),
  #            size = 5) +
  # geom_segment(data = LOH_type_xTx, 
  #            aes(x = low95CI_i, xend = up95CI_i, 
  #                y = y_pos_mean_i + (4 - as.numeric(Tx_ID)) * y_space_i, 
  #                yend = y_pos_mean_i + (4 - as.numeric(Tx_ID)) * y_space_i,
  #                color = Tx_ID),
  #            size = 1) +
  scale_fill_manual(values = txPal, guide = "none") +
  scale_color_manual(values = txPal, guide = "none") +
  scale_x_continuous(breaks = seq(0, max(iLOH_xn_count$n_iLOH + 1), 1)) +
  scale_y_continuous(breaks = seq(0, y_pos_mean_i, 5)) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_blank(),
        axis.title = element_blank())

y_pos_mean_t <- max(tLOH_xn_count$n)
y_space_t <- 2

tLOH_density <- tLOH_xn_count %>% 
  filter(Tx_ID != "D") %>%
  ggplot() +
  geom_col(aes(x = n_tLOH, y = n, fill = Tx_ID, color = Tx_ID), 
           alpha = 0.7,
           position = position_dodge(preserve = "single", width = 0.8), width = 0.8) +
  # geom_point(data = LOH_type_xTx, 
  #            aes(x = mean_tLOH, y = y_pos_mean_t + as.numeric(Tx_ID) * y_space_t,
  #                color = Tx_ID),
  #            size = 5) +
  # geom_segment(data = LOH_type_xTx, 
  #              aes(x = low95CI_t, xend = up95CI_t, 
  #                  y = y_pos_mean_t + as.numeric(Tx_ID) * y_space_t, 
  #                  yend = y_pos_mean_t + as.numeric(Tx_ID) * y_space_t,
  #                  color = Tx_ID),
  #              size = 1) +
  scale_fill_manual(values = txPal, guide = "none") +
  scale_color_manual(values = txPal, guide = "none") +
  scale_x_continuous(breaks = seq(0, max(iLOH_xn_count$n_iLOH + 1), 1)) +
  scale_y_continuous(breaks = seq(0, y_pos_mean_t, 5)) +
  theme(panel.grid = element_blank(),
        axis.text.y = element_blank(),
        axis.title = element_blank()) +
  coord_flip()


LOH_type_plot_1 <- iLOH_density + plot_spacer() + LOH_type_plot + tLOH_density + 
  plot_layout(ncol = 2, nrow = 2, 
    widths = c(4, 1), heights = c(1, 4)) 

ggsave(file.path(outIntDir, "iLOHvstLOH_counts_1_2022_05.png"), 
       plot = LOH_type_plot_1,
       device = "png",
       width = 11, height = 8.5, 
       units = "in",
       dpi = 600)

Length_data_in$err_type <- factor(ifelse(Length_data_in$is_error, 
                                       "single_error", 
                                       ifelse(Length_data_in$length == 1, 
                                              "single", "multiple")))

# Length_data_in$Tx_ID <- Recode_Tx_ID(Length_data_in$Tx_ID, "Tx_ID")

# Length_data_in %>%
#   group_by(Tx_ID) %>%
#   mutate(d = n)

###############################################################################
# Comparison of length distributions

LOH_length_ecd <- Length_data_in %>%
  ggplot() +
  # geom_polygon(data = ecd_sig_region, aes(x = POS/1000, y = n, fill = Tx_name)) +
  stat_ecdf(aes(x = log10(est_length), color = Tx_ID)) + 
  # geom_segment(data = ecd_bracket,
  #              aes(x = POS/1000, xend = POS/1000, y = y_lo, yend = y_hi)) +
  # geom_label(data = ad_label, aes(x = 10 , y = 0.85, label = sig, color = Tx_2), 
  #            size = 9, hjust = "left", label.size = 0, show.legend = F) +
  # geom_point(data = marker_set, aes(x = POS/1000, y = -0.05), shape = "|", size = 0.5) +
  # stat_ecdf(aes(x = POSi), color = "red", size = 0.3) + 
  scale_color_manual(values = txPal, name = "Strain") + 
  xlab("Chromosome Position (kbp)") + 
  ylab("Cumulative Fraction of LOH events") +
  scale_x_continuous(expand = c(0.01, 0.01)) +
  # facet_wrap(~rom_CHROM, ncol = 4, scales = "free_x") +
  theme(
    # panel.grid.minor.y = element_blank(), 
    # panel.grid.minor.x = element_blank(),
    # panel.background = element_rect(color="grey80"),
    text = element_text(size = 24))

LOH_length_ecd

LOH_length_type_wass <- wass_test(Length_data_in %>% filter(!isTerm) %>% pull(est_length),
                                  Length_data_in %>% filter(isTerm) %>% pull(est_length), 
                                  nboots = 5000)

LOH_length_type_wcx <- wilcox.test(Length_data_in %>% filter(!isTerm) %>% pull(est_length),
                                  Length_data_in %>% filter(isTerm) %>% pull(est_length))

all_LOH_length_wass <- data.frame(NULL)

for(ty in levels(Length_data_in$type)){
  type_LOH_length_wass <- data.frame(type = ty, Combo = with(Tx_ID_combo, paste(Tx_1, Tx_2, sep = "-")),
                                    Tx_1 = Tx_ID_combo[, 1], Tx_2 = Tx_ID_combo[, 2],
                                    W_stat = 0, W_pVal = 0, X_stat = 0, X_pVal = 0)
  n_c <- ncol(type_LOH_length_wass)
  for(tx in 1:nrow(Tx_ID_combo)){
    # tx = 1
    Tx_1 <- Tx_ID_combo$Tx_1[tx]
    Tx_2 <- Tx_ID_combo$Tx_2[tx]
    LOH_1 <- Length_data_in %>% 
      filter(Tx_ID == Tx_1, type == ty) %>% 
      select(Tx_ID, ID, est_length)
    LOH_2 <- Length_data_in %>% 
      filter(Tx_ID == Tx_2, type == ty) %>% 
      select(Tx_ID, ID, est_length)
    # AD_LOHcounts$Tx_ID <- droplevels(AD_LOHcounts$Tx_ID)
    length_wass <- wass_test(LOH_1$est_length, LOH_2$est_length, nboots = 5000)
    length_wcx <- wilcox.test(LOH_1$est_length, LOH_2$est_length)
    type_LOH_length_wass[tx, c(n_c - 3, n_c - 2)] <- length_wass[1:2]
    type_LOH_length_wass[tx, c(n_c - 1)] <- length_wcx$statistic
    type_LOH_length_wass[tx, c(n_c)] <- length_wcx$p.value
  }
  all_LOH_length_wass <- rbind(all_LOH_length_wass, type_LOH_length_wass)
}

all_LOH_length_wass
all_LOH_length_wass_slim <- all_LOH_length_wass %>% filter(Combo != "C-D")
all_LOH_length_wass_slim <- BHcorrection(all_LOH_length_wass_slim)
all_LOH_length_wass_slim$p_adjust <- p.adjust(all_LOH_length_wass_slim$pVal, method = "BH")

hist_group <- function(x, b = seq(1, 6.5, 0.2)) {
  return(data.frame(hist(log10(x), breaks = b, plot = F)[2:4]))
  # data.frame(breaks = h$breaks, mids = h$mids counts = h$counts, density = h$density)
}

length_hist <- Length_data_in %>% 
  group_by(type, Tx_ID) %>% 
  group_modify(~hist_group(.x$est_length), .keep = T) %>%
  mutate(freq = counts/sum(counts))

length_med_xTx <- Length_data_in %>% 
  group_by(type) %>% 
  # group_by(Tx_ID, type) %>% 
  summarize(m = median(est_length), .groups = "keep")

type_label <- data.frame(.x = 1.5, .y = c(0.21, 0.31),
                         type = factor(c("Interstitial", "Terminal")), 
                         l = factor(c("Interstitial LOH", "Terminal LOH")))

length_dist_hist <- #length_hist %>% 
  Length_data_in %>%
  # filter(Tx_name == "WT") %>%
  # filter(!is_error) %>%
  # filter(!isTerm) %>%
  # count(LOH_cmplx)
  # summarize(m = median(est_length))
  ggplot() + 
  geom_label(data = type_label, aes(x = .x, y = .y, label = l),
             label.size = 0, size = 9) +
  # geom_label(data = n_clones_LOH_xTx, aes(x = 6.25, y = 0.28, label = Tx_ID),
  #            label.size = 0, size = 9) +
  geom_histogram(aes(x = log10(est_length), y = stat(density*width),
                     fill = Tx_ID),
                 position = position_dodge2(width = 1, preserve = "single", padding = 0),
                 binwidth = 1/4, color = "grey20", size = 0.25) +
  geom_histogram(aes(x = log10(est_length), y = stat(density*width)),
                 binwidth = 1/4, color = "brown4", size = 0.75, fill = "white", alpha = 0) +
  # geom_vline(data = length_med_xTx, aes(xintercept = m), color = "blue3") +
  # geom_col(aes(x = mids, y = freq,
  #                    fill = Tx_ID), 
  #                position = position_dodge2(width = 1, preserve = "single", padding = 0),
  #                color = "grey20") + 
  scale_fill_manual(values = txPal, name = "Strain") +
  # scale_alpha_discrete(range = c(0.8, 0.3), name = "LOH Type") +
  scale_x_continuous(name = "LOH length (bp)", breaks = seq(2, 6, 1), 
                     labels = format_sci_10_log
                     ) +
  # guides(fill = "none") +
  ylab("Proportion of Events") +
  # facet_grid(Tx_ID~.) +
  facet_grid(type~., scales = "free_y") +
  # facet_grid(is_error~.) +
  # facet_grid(Tx_ID~type) +
  theme(legend.position = c(0.9, 0.92),
        text = element_text(size = 28),
        strip.text = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank())

length_dist_hist

ggsave(file.path(outIntDir, "LOH_lengthDist_hist_2022_05.png"), 
       plot = length_dist_hist,
       device = "png",
       width = 16, height = 12, 
       units = "in",
       dpi = 600)

Length_data_in %>% 
  filter(est_length > 10^6)
  filter(est_length > 10^4.1, est_length < 10^4.5, Tx_ID == "D") %>% 
  count(isTerm)

# Length_data_in$err_type <- factor(Length_data_in$err_type, 
#                                   levels = c("single_error", 
#                                              "single", "multiple"))

marker_type_dist <- Length_data_in %>% 
  count(Tx_name, err_type) %>% group_by(Tx_name) %>%
  mutate(fract_n = n/sum(n)) %>%
  ggplot() + 
  # geom_col(aes(x = Tx_name, y = fract_n, fill = err_type),
  #          color = "grey20") + 
  geom_col(aes(x = Tx_name, y = n, fill = err_type),
                      color = "grey20") + 
  scale_fill_manual(values = c( "brown3", "dodgerblue3", "grey40"), 
                    name = "markers") +
  xlab("Strain") +
  ylab("Number of LOH events") +
  theme(legend.position = "bottom",
        text = element_text(size = 14))

marker_type_dist

Length_data_in %>% count(Tx_ID) %>% pull(n) %>% min()
n_LOH_xTx <- all_LOHcounts_merge_NS %>% group_by(Tx_ID) %>% summarize(n = sum(n_LOH))
# n_error_xTx <- all_LOHcounts_merge %>% group_by(Tx_ID) %>% summarize(n = sum(n_LOH))
all_GTbounds_merge %>% 
  filter(GT != "0/1") %>%
  filter(length == 1) %>%
  count(Tx_ID, is_error) %>% group_by(Tx_ID) %>% summarize(percent_err = n/sum(n))
  pivot_wider(names_from = is_error, names_prefix = "err_", values_from = n) %>%
  mutate(percent_err = err_TRUE/(err_FALSE + err_TRUE))

  
  
lengthEC_dist_hist <- Length_data_in %>% 
  group_by(Tx_ID) %>%
  slice_sample(n = min(n_LOH_xTx$n), replace = F) %>%
  # filter(!is_error) %>% 
  # filter(!isTerm) %>%
  # group_by(isTerm) %>%
  # count(LOH_cmplx)
  # summarize(m = median(est_length))
  ggplot() + geom_histogram(aes(x = log10(est_length), fill = type),
                            binwidth = 0.25, color = "grey20") + 
  scale_fill_manual(values = c("grey40", "dodgerblue3","brown3"), 
                    name = "markers") +
  scale_x_continuous(name = "LOH length (bp)", breaks = seq(2, 6, 1), 
                     labels = formatC(10^seq(2, 6, 1), 1, format = "E")) +
  ylab("Number of LOH events") +
  # facet_grid(type~.) +
  facet_grid(Tx_ID~type) +
  theme(legend.position = "bottom",
        text = element_text(size = 14))

lengthEC_dist_hist

ggsave(file.path(outIntDir, "LOH_lengthDistEC_hist.png"), 
       plot = lengthEC_dist_hist,
       device = "png",
       width = 11, height = 8.5, 
       units = "in",
       dpi = 600)

lengthEC_dist_hist_xTx <- Length_data_in %>% 
  # filter(!is_error) %>%
  # group_by(isTerm) %>%
  # count(LOH_cmplx)
  # summarize(m = median(est_length))
  ggplot() + geom_histogram(aes(x = log10(est_length), fill = type),
                            binwidth = 0.25, color = "grey20") + 
  scale_fill_manual(values = c("grey30", "grey80"), 
                    name = "LOH Type") +
  scale_x_continuous(name = "LOH length (bp)", breaks = seq(2, 6, 1), 
                     labels = formatC(10^seq(2, 6, 1), 1, format = "E")) +
  ylab("Number of LOH events") +
  # facet_grid(type~.) +
  facet_grid(Tx_ID~.) +
  theme(legend.position = "bottom",
        text = element_text(size = 14))

lengthEC_dist_hist_xTx

ggsave(file.path(outIntDir, "lengthEC_dist_hist_xTx.png"), 
       plot = lengthEC_dist_hist_xTx,
       device = "png",
       width = 11, height = 8.5, 
       units = "in",
       dpi = 600)


###############################################################################
Length_data_cln <- Length_data_in # %>% 
  # filter(!is_error)
Length_data_N <- Length_data_cln %>% filter(Tx_name == "WT")
Length_data_H <- Length_data_cln %>% filter(Tx_name == "Cas9")
Length_data_F <- Length_data_cln %>% filter(Tx_name == "Drive")

i_N <- Length_data_N %>% filter(!isTerm) %>% pull(est_length)
i_H <- Length_data_H %>% filter(!isTerm) %>% pull(est_length)
i_F <- Length_data_F %>% filter(!isTerm) %>% pull(est_length)

t_N <- Length_data_N %>% filter(isTerm) %>% pull(est_length)
t_H <- Length_data_H %>% filter(isTerm) %>% pull(est_length)
t_F <- Length_data_F %>% filter(isTerm) %>% pull(est_length)

KS_NvH_i <- ks.test(i_N, i_H, alternative = "two.sided")
KS_NvH_i_p <- signif(KS_NvH_i$p.value, 3)

KS_NvF_i <- ks.test(i_N, i_F, alternative = "two.sided")
KS_NvF_i_p <- signif(KS_NvF_i$p.value, 3)

KS_HvF_i <- ks.test(i_H, i_F, alternative = "two.sided")
KS_HvF_i_p <- signif(KS_HvF_i$p.value, 3)

KS_NvH_t <- ks.test(t_N, t_H, alternative = "two.sided")
KS_NvH_t_p <- signif(KS_NvH_t$p.value, 3)

KS_NvF_t <- ks.test(t_N, t_F, alternative = "two.sided")
KS_NvF_t_p <- signif(KS_NvF_t$p.value, 3)

KS_HvF_t <- ks.test(t_H, t_F, alternative = "two.sided")
KS_HvF_t_p <- signif(KS_HvF_t$p.value, 3)

KS_label <- c(paste0("Kolmogorov-Smirnov Test",
                     "\niLOH",
                     "\n  WTvCas9 p = ", KS_NvH_i_p,
                     "\n  WTvsDrive p = ", KS_NvF_i_p,
                     "\n  Cas9vsDrive p = ", KS_HvF_i_p),
              paste0("\ntLOH",
                "\n  WTvCas9 p = ", KS_NvH_t_p,
                "\n  WTvsDrive p = ", KS_NvF_t_p,
                "\n  Cas9vsDrive p = ", KS_HvF_t_p))

# AD test #####################################################################
Length_NvH_i <- Length_data_cln %>% filter(Tx_name != "Drive", !isTerm)
Length_NvH_i$Tx_name <- factor(Length_NvH_i$Tx_name)
AD_NvH_i <- ad.test(est_length ~ Tx_name, data = Length_NvH_i)
AD_NvH_i_p <- signif(max(AD_NvH_i$ad[,3]))

Length_NvF_i <- Length_data_cln %>% filter(Tx_name != "Cas9", !isTerm)
Length_NvF_i$Tx_name <- factor(Length_NvF_i$Tx_name)
AD_NvF_i <- ad.test(est_length ~ Tx_name, data = Length_NvF_i)
AD_NvF_i_p <- signif(max(AD_NvF_i$ad[,3]), 3)

Length_HvF_i <- Length_data_cln %>% filter(Tx_name != "WT", !isTerm)
Length_HvF_i$Tx_name <- factor(Length_HvF_i$Tx_name)
AD_HvF_i <- ad.test(est_length ~ Tx_name, data = Length_HvF_i)
AD_HvF_i_p <- signif(max(AD_HvF_i$ad[,3]), 3)

Length_NvH_t <- Length_data_cln %>% filter(Tx_name != "Drive", isTerm)
Length_NvH_t$Tx_name <- factor(Length_NvH_t$Tx_name)
AD_NvH_t <- ad.test(est_length ~ Tx_name, data = Length_NvH_t)
AD_NvH_t_p <- signif(max(AD_NvH_t$ad[,3]))

Length_NvF_t <- Length_data_cln %>% filter(Tx_name != "Cas9", isTerm)
Length_NvF_t$Tx_name <- factor(Length_NvF_t$Tx_name)
AD_NvF_t <- ad.test(est_length ~ Tx_name, data = Length_NvF_t)
AD_NvF_t_p <- signif(max(AD_NvF_t$ad[,3]), 3)

Length_HvF_t <- Length_data_cln %>% filter(Tx_name != "WT", isTerm)
Length_HvF_t$Tx_name <- factor(Length_HvF_t$Tx_name)
AD_HvF_t <- ad.test(est_length ~ Tx_name, data = Length_HvF_t)
AD_HvF_t_p <- signif(max(AD_HvF_t$ad[,3]), 3)

AD_label <- c(paste0("Anderson-Darling Test",
                   "\n  WTvCas9 p = ", AD_NvH_i_p,
                   "\n  WTvsDrive p = ", AD_NvF_i_p,
                   "\n  Cas9vsDrive p = ", AD_HvF_i_p),
                   paste0(
                   "\n  WTvCas9 p = ", AD_NvH_t_p,
                   "\n  WTvsDrive p = ", AD_NvF_t_p,
                   "\n  Cas9vsDrive p = ", AD_HvF_t_p))

# AD test end #################################################################

all_med_len <- Length_data_cln %>% 
  group_by(type, Tx_ID) %>% 
  summarise(med_len = median(est_length))

label_df <- data.frame(est_length = c(100000, 50), cum_prop = 0.25, 
                       Tx_ID = "D", type = c("Interstitial", "Terminal"))

LOH_lengthDist <- Length_data_cln %>% 
  ggplot() + 
  geom_vline(data = all_med_len, aes(xintercept = med_len, color = Tx_ID), size = 0.25) +
  geom_label_repel(data = all_med_len, aes(x = med_len, y = 1, color = Tx_ID, label = round(med_len)), 
                   label.size = 0, box.padding = unit(0.3, "lines"), ylim = c(1.03, 1.07), key_glyph = "rect") +
  # geom_label(data = label_df, aes(x = est_length, y = cum_prop, label = KS_label),
  #            hjust = 0, label.size = 0) +
  stat_ecdf(aes(x = est_length, color = Tx_ID), pad = F) +
  scale_color_manual(values = txPal, name = "Drive type") +
  # scale_linetype_manual(name = "LOH type", values = c("solid", "dashed")) +
  scale_x_log10(breaks = c(1E0, 1E1, 1E2, 1E3, 1E4, 1E5, 1E6), limits = c(10, 2E6)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1.1)) +
  ylab("Cumulative proportion") + xlab("LOH Event Length (bp)") +
  theme(panel.grid.minor.x = element_blank(),
        panel.background = element_rect(color = "grey70"),
        text = element_text(size = 14),
        legend.position = "bottom") +
  facet_grid(type~., scales = "free_x")

LOH_lengthDist

ggsave(file.path(outIntDir, "LOH_lengthDist_plot_wide.png"), 
       plot = LOH_lengthDist,
       device = "png",
       width = 16, height = 9, 
       units = "in",
       dpi = 600)

Length_data_cln %>% 
  filter(CHROM == "05") %>%
  # filter(type == "Interstitial") %>%
  ggplot() + 
  geom_histogram(aes(x = log10(est_length), fill = Tx_ID), binwidth = 0.33) +
  scale_fill_manual(values = txPal, name = "Drive type") +
  # scale_y_continuous(expand = c(0, 0), limits = c(0, 1.1)) +
  # ylab("Cumulative proportion") + 
  xlab("Log(LOH Event Length (bp))") +
  theme(panel.grid.minor.x = element_blank(),
        panel.background = element_rect(color = "grey70"),
        text = element_text(size = 14),
        legend.position = "bottom") +
  facet_grid(Tx_ID~., scales = "free_x")


###############################################################################

# Counts of terminal and interstitial LOH events ------ #######################

LOH_type_counts <- all_LOHcounts_merge_EC %>% 
  group_by(Tx_name) %>% 
  summarise(sum_i = sum(n_Inter), sum_t = sum(n_Term), ratio = sum(n_Inter)/sum(n_Term))

# M <- as.table(rbind(c(762, 327, 468), c(484, 239, 477)))
# dimnames(M) <- list(gender = c("F", "M"),
#                     party = c("Democrat","Independent", "Republican"))
# Xsq <- chisq.test(M)

LOH_type_counts_table <- as.table(as.matrix(LOH_type_counts[, 2:3]))

chisq.test(LOH_type_counts_table[-2, ])

LOH_type_counts <- all_LOHcounts_EC %>% 
  select(Tx_name, ID, n_LOH, n_Term, n_Inter) %>% 
  pivot_longer(cols = c(n_Term, n_Inter), names_to = "type_name", values_to = "n")

LOH_type_counts$type <- ifelse(LOH_type_counts$type_name == "n_Term", "Terminal", "Interstitial")

LOH_type_counts <- as.data.frame(LOH_type_counts)

# Permutation test of the number of each LOH type
iLOH_counts_perm_N_F <- LOH_type_counts %>% 
  filter(type == "Interstitial") %>% 
  perm_test(., cat_var = "Tx_name", cat_names = c("WT", "Drive"), 
            response_var = "n", test_stat = "chisq_hist")
iLOH_counts_perm_N_H <- LOH_type_counts %>% 
  filter(type == "Interstitial") %>% 
  perm_test(., cat_var = "Tx_name", cat_names = c("WT", "Cas9"), 
            response_var = "n", test_stat = "chisq_hist")
iLOH_counts_perm_H_F <- LOH_type_counts %>% 
  filter(type == "Interstitial") %>% 
  perm_test(., cat_var = "Tx_name", cat_names = c("Cas9", "Drive"), 
            response_var = "n", test_stat = "chisq_hist")

iLOH_p <- round(c(iLOH_counts_perm_N_H$pVal, iLOH_counts_perm_N_F$pVal, 
                  iLOH_counts_perm_H_F$pVal), 3)
iLOH_is_sig <- ifelse(iLOH_p <= 0.05, "*", "")
iLOH_p_lab <- paste0("p = ", iLOH_p, iLOH_is_sig)

tLOH_counts_perm_N_F <- LOH_type_counts %>% 
  filter(type == "Terminal") %>% 
  perm_test(., cat_var = "Tx_name", cat_names = c("WT", "Drive"), 
            response_var = "n", test_stat = "chisq_hist")

tLOH_counts_perm_N_H <- LOH_type_counts %>% 
  filter(type == "Terminal") %>% 
  perm_test(., cat_var = "Tx_name", cat_names = c("WT", "Cas9"), 
            response_var = "n", test_stat = "chisq_hist")

tLOH_counts_perm_H_F <- LOH_type_counts %>% 
  filter(type == "Terminal") %>% 
  perm_test(., cat_var = "Tx_name", cat_names = c("Cas9", "Drive"), 
            response_var = "n", test_stat = "chisq_hist")

tLOH_p <- round(c(tLOH_counts_perm_N_F$pVal, tLOH_counts_perm_N_H$pVal, 
                  tLOH_counts_perm_H_F$pVal), 3)
tLOH_is_sig <- ifelse(tLOH_p <= 0.05, "*", "")
tLOH_p_lab <- paste0("p = ", tLOH_p, tLOH_is_sig)

iLOH_lab_y <- max(LOH_type_counts$n)

iLOHcount_xTx_boxplot <- LOH_type_counts %>% 
  ggplot() + 
  facet_grid(.~type) +
  geom_boxplot(aes(x = Tx_name, y = n, color = type), 
               width = 0.5, outlier.color = NA) +
  # geom_violin(aes(x = Tx_name, y = n, color = type), 
  #             width = 0.5) +
  geom_jitter(aes(x = Tx_name, y = n, color = type), alpha = 0.7, size = 2.5,
              position = position_jitterdodge(jitter.height = 0.1, 
                                              jitter.width = 0.3, 
                                              dodge.width = 0.5)) +
  # geom_bracket(xmin = c(1, 1, 2) - 0.15, xmax = c(2, 3, 3) - 0.15,
  #              y.position = c(iLOH_lab_y + 3, iLOH_lab_y + 6, iLOH_lab_y + 9),
  #              label = iLOH_p_lab,
  #              tip.length = 0.025, label.size = 4) +
  scale_color_manual(values = c("black", "dodgerblue3"), name = "LOH type") +
  ylim(c(NA, iLOH_lab_y + 10)) + ylab("Number of LOH events") + xlab("Drive type") +
  theme(text = element_text(size = 14))

iLOHcount_xTx_boxplot

ggsave(file.path(outIntDir, "LOHcount_TypexTx_boxplot_2021_12.png"), 
       plot = iLOHcount_xTx_boxplot,
       device = "png",
       width = 12, height = 10, 
       units = "in",
       dpi = 600)

LOH_type_count_means <- LOH_type_counts %>% 
  group_by(Tx_name, type) %>% summarize(m = median(n))
plot_bars <- data.frame(from = seq(0.5, 13.5, 2), to = seq(1.5, 14.5, 2))

LOH_type_counts %>% 
  ggplot() + 
  # geom_rect(data = plot_bars,
  #           aes(xmin = from, xmax = to,
  #               ymin = -Inf, ymax = Inf),
  #           alpha = 0.2, fill = "grey50") +
  geom_segment(data = LOH_type_count_means, 
               aes(x = m + 0.05*(as.numeric(Tx_name) - 2), xend = m + 0.05*(as.numeric(Tx_name) - 2), y = 0, yend = Inf, color = Tx_name),
             size= 0.4) +
  geom_histogram(aes(x = n, fill = Tx_name), 
                 # position = position_dodge2(preserve = "total", width = 1, padding = 0.2), 
                 binwidth = 1) +
  # geom_freqpoly(aes(x = n, color = Tx_name), binwidth = 1) +
  # geom_bracket(xmin = c(1, 1, 2) - 0.15, xmax = c(2, 3, 3) - 0.15,
  #              y.position = c(iLOH_lab_y + 3, iLOH_lab_y + 6, iLOH_lab_y + 9),
  #              label = iLOH_p_lab, 
  #              tip.length = 0.025, label.size = 4) +
  scale_color_manual(values = txPal, name = "Drive type") +
  scale_fill_manual(values = txPal, name = "Drive type") +
  scale_x_continuous(breaks = seq(0, max(LOH_type_counts$n), 1)) +
  # ylim(c(NA, iLOH_lab_y + 10)) + 
  xlab("# of LOH events") + ylab("Number of Clones") +
  facet_grid(Tx_name~type, scales = "free") + 
  theme(panel.grid.major.x = element_blank())


###############################################################################

LOH_type_counts_tx <- Length_data_in %>% 
  count(Tx_name, isTerm) %>% 
  pivot_wider(names_from = isTerm, names_prefix = "is_", values_from = n) %>%
  mutate(term = is_TRUE, inter = is_FALSE) %>% select(-c(is_TRUE, is_FALSE))

LOH_type_counts_tx$fTerm <- LOH_type_counts_tx$term/(LOH_type_counts_tx$term + LOH_type_counts_tx$inter)

# Distribution of number of occasions a site is converted among all clones ----
Length_data_in %>% 
  filter(length == 1) %>% 
  # select(start_POSi, GT) %>% 
  count(start_POSi, GT) %>% filter(n >= 2) %>% ggplot() + 
  geom_histogram(aes(x = n), binwidth = 1)

#####
# Length distribution by genotype ----
GT_med_len <- Length_data_in %>% 
  # filter(GT != "0/1") %>%
  group_by(isTerm, GT) %>% 
  summarise(med_len = round(median(est_length)))
  
LOH_lengthDist_GT <- Length_data_in %>% 
  # filter(GT != "0/1") %>%
  ggplot(aes(x = est_length, linetype = isTerm, color = GT)) + 
  geom_vline(data = GT_med_len, aes(xintercept = med_len, color = GT, linetype = isTerm), size = 0.25) +
  geom_label_repel(data = GT_med_len, aes(x = med_len, y = 1, color = GT, label = med_len), 
                   label.size = 0, box.padding = unit(0.5, "lines"), ylim = c(1, 1.05), key_glyph = "rect") +
  stat_ecdf(pad = F) +
  scale_color_manual(values = allelePal[c(1, 3)]) +
  scale_linetype_manual(name = "Type", 
                        values = c("solid", "dashed"), labels = c("iLOH", "tLOH")) +
  ylab("Cumulative proportion of LOH events") + xlab("LOH Event Length") + 
  ylim(c(0,1)) + scale_x_log10() +
  theme(legend.position = "bottom",
        text = element_text(size = 14))

LOH_lengthDist_GT

ggsave(file.path(outIntDir, "LOH_lengthDist_xGT.png"), 
       plot = LOH_lengthDist_GT,
       device = "png",
       width = 11, height = 8.5, 
       units = "in",
       dpi = 600)

length_med_xTx <- Length_data_in %>% 
  group_by(Tx_ID, isTerm) %>% summarize(m = round(median(est_length)))

LOH_lengthDist_log <- Length_data_in %>%
  ggplot(aes(x = est_length, color = Tx_ID, linetype = isTerm)) + 
  geom_vline(data = length_med_xTx, aes(xintercept = m, color = Tx_ID, linetype = isTerm), size = 0.25) +
  geom_label_repel(data = length_med_xTx, aes(x = m, y = 1, color = Tx_ID, label = m), 
                   label.size = 0, box.padding = unit(0.5, "lines"), ylim = c(1, 1.05), key_glyph = "rect") + 
  stat_ecdf(pad = F) +
  scale_color_manual(values = txPal, name = "Drive Type") +
  scale_linetype_manual(name = "Type", 
                        values = c("solid", "dashed"), labels = c("iLOH", "tLOH")) +
  ylab("Cumulative proportion of LOH events") + xlab("LOH Event Length (bp)") + 
  scale_x_log10() +
  theme(legend.position = "bottom",
        text = element_text(size = 14))

LOH_lengthDist_log

ggsave(file.path(outIntDir, "LOH_lengthDist_log.png"), 
       plot = LOH_lengthDist_log,
       device = "png",
       width = 11, height = 8.5, 
       units = "in",
       dpi = 600)

#####


# Average LOH length vs no of LOH for each clone --------------

LOH_medianLength <- Length_data_cln %>% 
  group_by(Tx_ID, ID) %>% 
  summarise(med_len = median(est_length),
            mean_len = mean(est_length),
            med_log_len = mean(log10(est_length)))

LOH_n <- Length_data_cln %>%
  # group_by(Tx_ID, ID) %>% 
  count(Tx_ID, ID)

LOH_medianLength <- merge(LOH_medianLength, LOH_n, by = c("Tx_ID", "ID"), all = T)

LOH_len_n_lm <- lm(n ~ log10(med_len), data = LOH_medianLength)
summary(LOH_len_n_lm)

LOH_lenVsN <- ggplot(LOH_medianLength, aes(log10(med_len), n)) + 
  # geom_hline(aes(yintercept = 0), size = 0.15) +
  geom_smooth(method = "lm", fill = "grey80") +
  geom_jitter(aes(color = Tx_ID), height = 0.2, width = 0) +
  xlab("log(median LOH length)") + ylab("Number of LOH Events") +
  xlim(c(1, NA)) + #scale_x_log10() +
  scale_color_manual(values = txPal)

LOH_lenVsN

ggsave(file.path(outIntDir, "LOH_lenVsN.png"), 
       plot = LOH_lenVsN,
       device = "png",
       width = 10, height = 10, 
       units = "in",
       dpi = 300)


# UNUSED ######################################################################

# For each LOH type, iLOH and tLOH, distribution of positions along chromosome ----

evo_homRuns_pos <- all_homRuns %>% filter(Rep != "00") %>% 
  select(Tx_name, ID, GT, isTerm, chrom_n, startPOSi, endPOSi)

evo_homRuns_pos$startPOS <- 0
evo_homRuns_pos$endPOS <- 0

for (ch in 1:nrow(chrom_bound_BY)) {
  # ch = 1
  idx <- evo_homRuns_pos$chrom_n == ch
  evo_homRuns_pos$startPOS[idx] <- evo_homRuns_pos$startPOSi[idx] - chrom_bound_BY[ch,1] + 1
  evo_homRuns_pos$endPOS[idx] <- evo_homRuns_pos$endPOSi[idx] - chrom_bound_BY[ch,1] + 1
}


# T-------L----O-----------T

evo_homRuns_pos$dCent <- 0


for (ch in 1:nrow(centrom_df)) {
  # ch = 1
  idx <- evo_homRuns_pos$chrom_n == ch & evo_homRuns_pos$isTerm == F
  dwnStrmIdx <- evo_homRuns_pos$startPOS[idx] > centrom_df$POS[ch]
  upStrmIdx <- evo_homRuns_pos$startPOS[idx] < centrom_df$POS[ch] & evo_homRuns_pos$endPOS[idx] < centrom_df$POS[ch]
  crsStrmIdx <- evo_homRuns_pos$startPOS[idx] < centrom_df$POS[ch] & evo_homRuns_pos$endPOS[idx] > centrom_df$POS[ch]
  df <- evo_homRuns_pos[idx,]
  df[crsStrmIdx,]
  evo_homRuns_pos$dCent[idx] <- 
}

# iLOH that cross the centromere 
# chrom 4(7), 6(1), 14(1), 

# Plot location and number of LOHs across genome ---------

# all_GTruns_cln

all_LOHbounds$midPOSi <- (all_LOHbounds$est_start + all_LOHbounds$est_end)/2

all_LOHbounds %>% 
  filter(GT != "het") %>%
  ggplot() +
  scale_x_continuous(labels=as.character(chrom_bound_BY$CHROM),
                     breaks=(chrom_bound_BY$Start+chrom_bound_BY$End)/2) +
  annotate(geom="rect", xmin=chrom_bound_BY$Start[c(TRUE, FALSE)],
           xmax=chrom_bound_BY$End[c(TRUE, FALSE)], 
           ymin=-Inf, ymax=Inf, alpha=0.05, color = "grey90") +
  geom_point(aes(x = midPOSi, y = ID, color= GT), size = 0.9, shape = 15) +
  geom_segment(aes(x = est_start, xend = est_end, 
                   y = ID, yend = ID, color = GT), size = 1) +
  scale_color_manual(values = allelePal[c(1, 3)]) +
  theme(panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(color = "grey50", size = 0.15)) +
  facet_wrap(.~Tx_name, ncol = 1, scales = "free")

all_ancHet %>% filter(POSi >= 5828400, POSi < 5829000) %>%
  ggplot() +
  geom_point(aes(x = POSi, y = ID, color= GT), size = 0.9, shape = 15) +
  scale_color_manual(values = allelePal) +
  facet_wrap(.~Tx_name, ncol = 1, scales = "free")

all_ancHet %>% filter(POSi >= 5828400, POSi < 5829000) %>%
  ggplot() +
  geom_abline() +
  geom_point(aes(x = BY_DP_BYcall, y = RM_DP_BYcall, color= GT_BYcall)) +
  scale_color_manual(values = allelePal)

all_ancHet %>% filter(POSi >= 5828400, POSi < 5829000) %>%
  ggplot() +
  geom_abline() +
  geom_point(aes(x = BY_DP_RMcall, y = RM_DP_RMcall, color= GT_RMcall)) +
  scale_color_manual(values = allelePal)

# Plot LOH events ordered by position along chromosomes ------
all_LOHbounds_order <- all_LOHbounds %>% 
  filter(GT != "het") %>%
  group_by(Tx_name, CHROM) %>% 
  arrange(est_start) %>% 
  mutate(order = row_number())

LOH_stack_plot <- all_LOHbounds_order %>%
  ggplot() +
  scale_x_continuous(labels=as.character(chrom_bound_BY$CHROM),
                     breaks=(chrom_bound_BY$Start+chrom_bound_BY$End)/2) +
  annotate(geom="rect", xmin=chrom_bound_BY$Start[c(TRUE, FALSE)],
           xmax=chrom_bound_BY$End[c(TRUE, FALSE)], 
           ymin=-Inf, ymax=Inf, alpha=0.05, color = "grey90") +
  geom_segment(aes(x = est_start, xend = est_end, 
                   y = order, yend = order, color = GT), 
               alpha = 0.4, size = 1) +
  geom_segment(aes(x = midPOSi - 5000, xend = midPOSi + 5000, 
                   y = order, yend = order, color = GT), 
               alpha = 1, size = 1) +
  geom_point(data = centrom_df, aes(x = POSi, y = -5e-6), color = "blue4", shape = 17, size = 1.5) +
  # geom_point(aes(x = midPOSi, y = order, color= GT), size = 1.5, shape = "|") +
  scale_color_manual(values = allelePal[c(1, 3)]) +
  theme(panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(color = "grey50", size = 0.15)) +
  facet_wrap(.~Tx_name, ncol = 1, scales = "free")

LOH_stack_plot

ggsave(file.path(outIntDir, "LOH_stack_plot.png"), 
       plot = LOH_stack_plot,
       device = "png",
       width = 16, height = 9, 
       units = "in",
       dpi = 300)

all_LOHbounds_order %>%
  ggplot() +
  scale_x_continuous(labels=as.character(chrom_bound_BY$CHROM),
                     breaks=(chrom_bound_BY$Start+chrom_bound_BY$End)/2) +
  annotate(geom="rect", xmin=chrom_bound_BY$Start[c(TRUE, FALSE)],
           xmax=chrom_bound_BY$End[c(TRUE, FALSE)], 
           ymin=-Inf, ymax=Inf, alpha=0.05, color = "grey90") +
  geom_segment(aes(x = est_start, xend = est_end, 
                   y = ID, yend = ID, color = GT), 
               alpha = 0.7, size = 1) +
  geom_segment(aes(x = midPOSi - 5000, xend = midPOSi + 5000, 
                   y = ID, yend = ID, color = GT), 
               alpha = 1, size = 1) +
  geom_point(data = centrom_df, aes(x = POSi, y = -5e-6), color = "blue4", shape = 17, size = 1.5) +
  # geom_point(aes(x = midPOSi, y = order, color= GT), size = 1.5, shape = "|") +
  scale_color_manual(values = allelePal[c(1, 3)]) +
  theme(panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(color = "grey50", size = 0.15)) +
  facet_wrap(.~Tx_name, ncol = 1, scales = "free")


all_ancHet %>% 
  filter(CHROM == "10") %>% 
  ggplot(aes(x = POSi)) + 
  geom_histogram(aes(fill = Tx_name), 
                 binwidth = 10000,
                 position = position_dodge2(preserve = "total")) +
  scale_fill_manual(values = txPal)


evo_homRuns_srt <- evo_homRuns_cln %>% dplyr::arrange(Tx, chrom_n, midPOSi), GT)

win <- 50000
chrom_wnds <- c()
for(ch in 1:16) {
  chrom_wnds <- c(chrom_wnds, seq(from = chrom_bound_BY[ch, 1], to = chrom_bound_BY[ch, 2], by = win), chrom_bound_BY[ch, 2])
}

evo_homRuns_out <- data.frame(NULL)
evo_homRuns_srt$y_coord <- 1
evo_homRuns_srt$wndw <- 0
for(tx in levels(evo_homRuns_srt$Tx)) {
  # tx = "N"
  tx_homRuns <- evo_homRuns_srt[evo_homRuns_srt$Tx == tx,]
  for(w in 1:(length(chrom_wnds)-1)) {
    # w = 31
    widx <- tx_homRuns$midPOSi >= chrom_wnds[w] & tx_homRuns$midPOSi <= chrom_wnds[w+1]
    i_seq <- 1:sum(widx)
    tx_homRuns[widx,] <- tx_homRuns[widx,] %>% arrange(GT)
    tx_homRuns$y_coord[widx] <- i_seq
    tx_homRuns$wndw[widx] <- (chrom_wnds[w] + chrom_wnds[w+1] - 2)/2
  }
  evo_homRuns_out <- rbind(evo_homRuns_out, tx_homRuns)
}


evo_homRuns_out <- data.frame(NULL)
evo_homRuns_srt$y_coord <- 1
evo_homRuns_srt$done <- 0
i = 1
for(tx in levels(evo_homRuns_srt$Tx)) {
  # tx = "N"
  tx_homRuns <- evo_homRuns_srt[evo_homRuns_srt$Tx == tx,]
  tx_df <- data.frame(NULL)
  for(ch in 1:16) {
    chr_homRuns <- tx_homRuns[tx_homRuns$chrom_n == ch,]
    for (l in 1:(nrow(chr_homRuns)-1)) {
      # l = 1
      if (chr_homRuns$done[l] == 0) {
        clsIdx <- abs(chr_homRuns$midPOSi[l] - chr_homRuns$midPOSi) < win
        i_seq <- 1:sum(clsIdx)
        chr_homRuns$y_coord[clsIdx] <- i_seq
        chr_homRuns$done[clsIdx] <- 1
      }
    }
    tx_df <- rbind(tx_df, chr_homRuns)
  }
  
  evo_homRuns_out <- rbind(evo_homRuns_out, tx_df)
}

evo_homRuns_out$startPOS <- 0
evo_homRuns_out$midPOS <- 0
evo_homRuns_out$endPOS <- 0
evo_homRuns_out$winPOS <- 0


for (ch in 1:nrow(chrom_bound_BY)) {
  # ch = 1
  idx <- evo_homRuns_out$chrom_n == ch
  evo_homRuns_out$startPOS[idx] <- evo_homRuns_out$startPOSi[idx] - chrom_bound_BY[ch,1] + 1
  evo_homRuns_out$midPOS[idx] <- evo_homRuns_out$midPOSi[idx] - chrom_bound_BY[ch,1] + 1
  evo_homRuns_out$endPOS[idx] <- evo_homRuns_out$endPOSi[idx] - chrom_bound_BY[ch,1] + 1
  evo_homRuns_out$winPOS[idx] <- evo_homRuns_out$wndw[idx] - chrom_bound_BY[ch,1] + 1
}

evo_homRuns_out$Tx_name <- Recode(evo_homRuns_out$Tx, "WT" = "N", "Cas9" = "H", "Drive" = "F")
evo_homRuns_out$Tx_name <- factor(evo_homRuns_out$Tx_name, levels = c("WT", "Cas9", "Drive"))

# This is close to the plot that we want, but the x scale cannot be adjusted
evo_homRuns_out %>%
  ggplot() + 
  geom_dotplot(aes(x = midPOS), binwidth = 50000) + facet_grid(Tx~chrom_n)
  
  
nLOHxPOS_SW <- evo_homRuns_out %>% 
  ggplot() + 
  # scale_x_continuous(labels=as.character(chrom_bound_BY$CHROM),
  #                    breaks=(chrom_bound_BY$Start+chrom_bound_BY$End)/2) +
  xlab("Chromosome") + ylab("# LOH events") + 
  theme(panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank(),
        axis.text.x = element_blank(),
        panel.background = element_rect(color="grey80"),
        strip.text.y = element_text(size = 11),
        panel.spacing = unit(0.1, "lines")) +
  # annotate(geom="rect", xmin=chrom_bound_BY$Start[c(TRUE, FALSE)],
  #          xmax=chrom_bound_BY$End[c(TRUE, FALSE)], ymin=-Inf, ymax=Inf, alpha=0.05, color = "grey90") +
  geom_point(aes(x = winPOS, y = y_coord, color = GT), size = 0.5) + 
  scale_color_manual(values = c("orange1", "blue2"), labels = c("BY", "RM")) +
  #geom_segment(aes(x = startPOSi, xend = endPOSi, y = y_coord, yend = y_coord), size = 0.5) + 
  facet_grid(Tx_name~chrom_n, space = "free_x", scales = "free_x", switch = "x")

ggsave(file.path(outIntDir, "nLOHxPOS_SW.png"), 
       plot = nLOHxPOS_SW,
       device = "png",
       width = 16, height = 9, 
       units = "in",
       dpi = 300)


WT_nLOH_SW <- SliderCalc(evo_homRuns_out[evo_homRuns_out$Tx == "N",], data_col = "midPOSi", index_col = "midPOSi", 
           factor_col = "chrom_n", 
           window_size = 100000, 
                       slide_interval, summary_stat = length)

WT_nLOH_SW %>%
  ggplot() + geom_col(aes(x = start, y = n_elements)) +
  scale_x_continuous(labels=as.character(chrom_bound_BY$CHROM),
                     breaks=(chrom_bound_BY$Start+chrom_bound_BY$End)/2) +
  annotate(geom="rect", xmin=chrom_bound_BY$Start[c(TRUE, FALSE)],
           xmax=chrom_bound_BY$End[c(TRUE, FALSE)], ymin=-Inf, ymax=Inf, alpha=0.2) +
  theme(panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank(),
        panel.background = element_rect(color="grey80"),
        strip.text.y = element_text(size = 11))



