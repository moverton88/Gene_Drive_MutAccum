
# Genomic iLOH and tLOH rates among strains


###############################################################################
# Proportion of iLOHs vs tLOHs different among pairs of strains ###############
# requires LOH_type_xTx from 07_LOH_length_analysis

LOHcounts_in <- all_LOHcounts_merge_NS


W_C_inter_perm <- LOHcounts_in %>% 
  # filter(!Line %in% c("H_F", "H_H")) %>%
  perm_test(df=., cat_var ="Tx_ID", cat_names = c("W", "C"), response_var = "n_Inter", 
            n_perms = 10000, alpha = 0.01, alt_hyp = "two-tailed", test_stat = mean, rtrn = "all") %>%
  as.data.frame() %>% mutate(Tx_ID = "C", type = "Interstitial")

W_D_inter_perm <- LOHcounts_in %>% 
  perm_test(df=., cat_var ="Tx_ID", cat_names = c("W", "D"), response_var = "n_Inter", 
            n_perms = 10000, alpha = 0.01, alt_hyp = "two-tailed", test_stat = mean, rtrn = "all") %>%
  as.data.frame() %>% mutate(Tx_ID = "D", type = "Interstitial")

W_C_term_perm <- LOHcounts_in %>% 
  # filter(!Line %in% c("H_F", "H_H")) %>%
  perm_test(df=., cat_var ="Tx_ID", cat_names = c("W", "C"), response_var = "n_Term", 
            n_perms = 10000, alpha = 0.01, alt_hyp = "two-tailed", test_stat = mean, rtrn = "all") %>%
  as.data.frame() %>% mutate(Tx_ID = "C", type = "Terminal")

W_D_term_perm <- LOHcounts_in %>% 
  perm_test(df=., cat_var ="Tx_ID", cat_names = c("W", "D"), response_var = "n_Term", 
            n_perms = 10000, alpha = 0.01, alt_hyp = "two-tailed", test_stat = mean, rtrn = "all") %>%
  as.data.frame() %>% mutate(Tx_ID = "D", type = "Terminal")

type_perm <- rbind(W_C_inter_perm, W_D_inter_perm, W_C_term_perm, W_D_term_perm)
type_perm$pVal %>% p.adjust()
type_sig <- type_perm %>% filter(pVal < 0.05)

type_z_prop <- LOH_type_xTx %>% 
  select(Tx_ID, n_iLOH, n_tLOH) %>% 
  proportion_z_score(grp = "Tx_ID") %>%
  filter(Tx_ID_null == "W")

fisher.test(t(LOH_type_xTx[c(3, 1), 2:3]))
fisher.test(t(LOH_type_xTx[c(3, 2), 2:3]))

tLOH_Tx_type_counts <- all_LOHcounts_merge_NS %>% 
  count(Tx_ID, n_Term) %>% rename(n_LOH = n_Term)

tLOH_Tx_type_counts <- fill_zeros(tLOH_Tx_type_counts, 
                                       group_col = "Tx_ID", 
                                       factor_col = "n_LOH", 
                                       fill_col = "n")
tLOH_Tx_type_counts$type <- "Terminal"


iLOH_Tx_type_counts <- all_LOHcounts_merge_NS %>% 
  count(Tx_ID, n_Inter) %>% rename(n_LOH = n_Inter)
iLOH_Tx_type_counts <- fill_zeros(iLOH_Tx_type_counts, 
                             group_col = "Tx_ID", 
                             factor_col = "n_LOH", 
                             fill_col = "n")
iLOH_Tx_type_counts$type <- "Interstitial"
LOH_Tx_type_counts <- rbind(iLOH_Tx_type_counts, tLOH_Tx_type_counts)
LOH_Tx_type_counts$type <- factor(LOH_Tx_type_counts$type)


# Plot center-out histogram of iLOH and tLOH count !!!!!!!####################

n_muts_lim <- max(LOH_Tx_type_counts$n) + 1
y_gap <- 0.05

LOH_type_means <- all_LOHcounts_merge_NS %>% 
  group_by(Tx_ID) %>% 
  summarize(Interstitial = mean(n_Inter),
            Terminal = mean(n_Term)) %>%
  pivot_longer(cols = -Tx_ID, names_to = "type", values_to = "m")

LOH_Tx_type_counts$x_left <- -LOH_Tx_type_counts$n/2
LOH_Tx_type_counts$x_right <- LOH_Tx_type_counts$n/2
LOH_Tx_type_counts$y_bot <- LOH_Tx_type_counts$n_LOH - 1/2 + y_gap
LOH_Tx_type_counts$y_top <- LOH_Tx_type_counts$n_LOH + 1/2 - y_gap

LOH_Tx_type_counts$text_pos <- ifelse(LOH_Tx_type_counts$n == 0 | LOH_Tx_type_counts$n > 1, 
                                 0, 
                                 LOH_Tx_type_counts$x_right + 1)

LOH_Tx_type_counts$text_color <- ifelse(LOH_Tx_type_counts$n > 1, 
                                              "white",
                                              txPal[as.numeric(LOH_Tx_type_counts$Tx_ID)])

text_backing <- LOH_Tx_type_counts %>% 
  filter(text_color != "white") %>% 
  select(Tx_ID, type, n_LOH, text_pos)
text_backing$x_min <- text_backing$text_pos - 1.75
text_backing$x_max <- text_backing$text_pos + 1.75
text_backing$y_min <- text_backing$n_LOH - 0.3
text_backing$y_max <- text_backing$n_LOH + 0.3

LOH_type_counts_rect_plot <- LOH_Tx_type_counts %>%
  ggplot() + 
  # geom_hline(data = SNMs_final_counts_stats, aes(yintercept = mean_n, color = Tx_name), size = 0.75) +
  geom_rect(data = text_backing, 
            aes(xmin = x_min, xmax = x_max, 
                ymin = y_min, ymax = y_max), fill = "white") +
  geom_rect(aes(xmin = x_left, xmax = x_right, ymin = y_bot, ymax = y_top, fill = Tx_ID)) +
  geom_crossbar(data = LOH_type_means, 
                aes(x = 0, y = m, ymin = m, ymax = m), 
                size = 0.1, width = n_muts_lim) +
  geom_text(aes(x = text_pos, y = n_LOH, 
                label = n), color = LOH_Tx_type_counts$text_color, size = 7) +
  scale_x_continuous(breaks = 0, label = NULL, limits = c(-n_muts_lim/2 - 2, n_muts_lim/2 + 2),
                     name = "Number of Clones") +
  scale_y_continuous(breaks = seq.int(0, max(LOH_Tx_type_counts$n, 1)), 
                     name = "Number of LOH events") +
  scale_fill_manual(values = txPal) + 
  scale_color_manual(values = txPal) +
  # coord_flip() +
  theme(legend.position = "none",
        panel.grid.minor.y = element_blank(),
        plot.margin = unit(c(t = 5, r = 5, b = 5, l = 8), "mm"),
        axis.title.x = element_text(vjust = 0),
        # axis.title.y = element_text(vjust = 3),
        strip.text.y = element_text(size = 22),
        strip.text.x = element_text(size = 28),
        text = element_text(size = 28)) +
  facet_grid(type~Tx_ID)

LOH_type_counts_rect_plot

ggsave(file.path(outIntDir, "LOH_type_counts_rect_plot_2022_04.png"),
       plot = LOH_type_counts_rect_plot,
       device = "png",
       width = 11, height = 8.5,
       units = "in",
       dpi = 600)

###############################################################################
# /bp rate of iLOH and tLOH given chromosome length ####################
LOHcounts_ID_CHROM_type <- POSi_data_in %>% 
  count(type, Tx_ID, ID, CHROM, .drop = F) %>%
  rename(n_LOH = n)

LOHcounts_CHROM_type_mean <- POSi_data_in %>% 
  count(type, Tx_ID, CHROM, .drop = F)  %>%
  rename(total_LOH = n)

LOHcounts_CHROM_type_mean$chrom_length <- rep(chrom_lengths_BY, 6)

LOH_CHROM_type_CI_list <- LOHcounts_ID_CHROM_type %>% 
  group_by(type, Tx_ID, CHROM) %>%
  group_modify(~ data.frame(boot.ci(boot(.x$n_LOH,
                                         statistic = sum.fun, R = 1000), 
                                    type = "perc")$percent[4:5] %>% t())) %>%
  rename(CI_95lo = X1, CI_95up = X2) %>% 
  select(CI_95lo, CI_95up)

LOHcounts_CHROM_type_mean <- merge(LOHcounts_CHROM_type_mean, 
                                   LOH_CHROM_type_CI_list,
                                   by = c("type", "Tx_ID", "CHROM"))

LOHcounts_CHROM_type_mean <- LOHcounts_CHROM_type_mean %>% 
  mutate(CHROM_bp_rate = total_LOH/chrom_length,
         CHROM_rate_CI_lo = CI_95lo/chrom_length, 
         CHROM_rate_CI_up = CI_95up/chrom_length)

LOHcounts_CHROM_type_mean <- LOHcounts_CHROM_type_mean %>% 
  mutate(bp_rate = CHROM_bp_rate / sum(n_clones_LOH_xTx$n) / n_gens,
         rate_CI_lo = CHROM_rate_CI_lo / sum(n_clones_LOH_xTx$n) / n_gens, 
         rate_CI_up = CHROM_rate_CI_up / sum(n_clones_LOH_xTx$n) / n_gens)

# LOHcounts_CHROM_type_mean$bp_rate <- operate_by_factor_match(
#   n_clones_LOH_xTx,
#   LOHcounts_CHROM_type_mean[, c("Tx_ID", "CHROM_bp_rate")],
#   .fun = function(x, y) y/x/n_gens)
# LOHcounts_CHROM_type_mean$bp_CI_lo <- operate_by_factor_match(
#   n_clones_LOH_xTx,
#   LOHcounts_CHROM_type_mean[, c("Tx_ID", "CHROM_rate_CI_lo")],
#   .fun = function(x, y) y/x/n_gens)
# LOHcounts_CHROM_type_mean$bp_CI_up <- operate_by_factor_match(
#   n_clones_LOH_xTx,
#   LOHcounts_CHROM_type_mean[, c("Tx_ID", "CHROM_rate_CI_up")],
#   .fun = function(x, y) y/x/n_gens)

iLOHcounts_CHROM_lm <- lm(bp_rate ~ chrom_length, 
                          data = LOHcounts_CHROM_type_mean %>% filter(type == "Interstitial"))
summary(iLOHcounts_CHROM_lm)

tLOHcounts_CHROM_lm <- lm(bp_rate ~ chrom_length, 
                          data = LOHcounts_CHROM_type_mean %>% filter(type == "Terminal"))
summary(tLOHcounts_CHROM_lm)

LOHcounts_CHROM_type_mean %>% 
  ggplot() + geom_point(aes(x = chrom_length, y = bp_rate, color = Tx_ID)) + facet_grid(type~.)

LOHcounts_CHROM_type_mean %>% 
  group_by(type) %>% 
  summarize(range_LOH = max(bp_rate)/min(bp_rate))

lm_coeff_list <- LOHcounts_CHROM_type_mean %>% 
  group_by(type) %>% 
  group_map(~ lm(bp_rate ~ chrom_length, data = .x),
            .keep = T)

lm_coeff_coeff <- lapply(lm_coeff_list, function(x) c(b = unname(x$coefficients[1]), 
                                                      m = unname(x$coefficients[2]),
                                                      rsq = summary(x)$adj.r.squared,
                                                      p = summary(x)$coefficients[2, 4]))

get_rsq_labels <- function(x) {
  # i_in <- format_sci_10(x[1])
  # s_in <- format_sci_10(x[2])
  # r_in <- format(x[3], digits = 2)
  i_in <- format(x[1], digits = 2)
  s_in <- format(x[2], digits = 2)
  r_in <- format(x[3], digits = 2)
  as.character(
    as.expression(
      substitute(y == i + s %*% x* "\n" ~~R^2~"="~r2,
                 list(i = i_in,
                      s = s_in,
                      r2 = r_in))
    )
  )
}


lm_labels <- lapply(lm_coeff_coeff, function(x) get_rsq_labels(x))

lm_coeff_df <- do.call(rbind, lm_coeff_coeff) %>% as.data.frame()
lm_labels_df <- do.call(rbind, lm_labels) %>% as.data.frame()
lm_coeff_df <- cbind(lm_coeff_df, label = lm_labels_df[, 1])
lm_coeff_df$type <- factor(c("Interstitial", "Terminal"))
# lm_coeff_df$Tx_ID <- factor(rep(Tx_ID_levels, 2), levels = Tx_ID_levels)
colnames(lm_coeff_df)[1:2] <- c("b", "m")

lm_coeff_df <- lm_coeff_df %>% 
  mutate(.x = min(chrom_lengths_BY), 
         .xend = max(chrom_lengths_BY), 
         .y = min(chrom_lengths_BY) * m + b, 
         .yend = max(chrom_lengths_BY) * m + b,
         l.x = max(chrom_lengths_BY) * 0.8,
         l.y = (max(.y) + 3E-10),
         eq = paste0("y == ", format_sci_10(b), " ", format_sci_10(m), "%.% x"),
         rsq_l = paste0(" \nR^2 == ", round(rsq, 3)),
         p_l = paste0(" \nP == ", round(p, 3))
         # lm_label = paste0(eq, "\nAdj R^2 = ", round(rsq, 3),
         #        "\nP = ", round(p, 3)),
         )
lm_coeff_df <- lm_coeff_df %>% mutate(lbl = paste0(eq, rsq_l, p_l))
# lm_label <- paste0(eq, "\nAdj R^2 = ", round(lm_coeff_df$rsq, 3), 
#                    "\nP = ", round(lm_coeff_df$p, 3))
# 
# lm_coeff_df$lm_label <- lm_label


BP_rate_chrom_length <- LOHcounts_CHROM_type_mean %>% 
  group_by(CHROM) %>%
  ggplot() +
  geom_point(aes(x = chrom_length, y = bp_rate), size = 3) +
  geom_label(data = lm_coeff_df, aes(x = l.x, y = l.y, label = eq),
             size = 6, hjust = 0, label.size = 0, show.legend = F, parse = T) +
  geom_label(data = lm_coeff_df, aes(x = l.x, y = l.y * 0.85, label = rsq_l),
             size = 6, hjust = 0, label.size = 0, show.legend = F, parse = T) +
  geom_label(data = lm_coeff_df, aes(x = l.x, y = l.y * 0.7, label = p_l),
             size = 6, hjust = 0, label.size = 0, show.legend = F, parse = T) +
  # annotate(geom = "label", x = lm_coeff_df$.x, y = lm_coeff_df$.y, label = lm_coeff_df$lm_lb, parse = F) +
  geom_segment(data = lm_coeff_df, aes(x = .x, xend = .xend, y = .y, yend = .yend), size = 1) +
  scale_color_manual(values = txPal, name = "Strain") +
  scale_x_continuous(labels = format_sci_10,
                     limits = c(0, NA), name = "Chromosome length (bp)") +
  scale_y_continuous(labels = format_sci_10,
                     limits = c(0, NA), name = "LOH event rate (/bp/gen)") +
  facet_wrap(~type, ncol = 1) +
  theme(text = element_text(size = 28),
        legend.position = c(0.93, 0.93))

BP_rate_chrom_length

ggsave(file.path(outIntDir, "BP_rate_chrom_length_2022_04.png"), 
       plot = BP_rate_chrom_length,
       device = "png",
       width = 16, height = 9, 
       units = "in",
       dpi = 600)


# Expected chromosome rate given genomic rate and chromosome size
expected_CHROM_rate <- lapply(chrom_lengths_BY_df$chrom_length,
                              function(y) sum(mean_LOHrate$sum_LOH)*(y/g_length)) %>%
  do.call("rbind", .) %>%
  as.data.frame() %>% rename(n_LOH = V1)

expected_CHROM_type_rate <- LOHcounts_CHROM_type_mean %>% 
  select(Tx_ID, CHROM, type, total_LOH) %>% 
  pivot_wider(names_from = "type", values_from = "total_LOH") %>%
  rename(total_iLOH = Interstitial, total_tLOH = Terminal) %>%
  arrange(Tx_ID, CHROM)

# f_iLOH <- tx_LOHrate %>% ungroup() %>% summarize(f = sum(total_iLOH)/sum(total_LOH))
expected_CHROM_type_rate <- lapply(tx_LOHrate$total_iLOH, 
                                function(x) sapply(chrom_lengths_BY_df$chrom_length,
                                                   function(y) x * (y/g_length))) %>% 
  as.data.frame(col.names = Tx_ID_levels) %>% 
  mutate(CHROM = chrom_lengths_BY_df$CHROM) %>% 
  pivot_longer(cols = 1:3, names_to = "Tx_ID", values_to = "expected_iLOH") %>%
  mutate(Tx_ID = factor(Tx_ID, Tx_ID_levels)) %>%
  arrange(Tx_ID, CHROM) %>% 
  select(expected_iLOH) %>%
  cbind(expected_CHROM_type_rate, .)

expected_CHROM_type_rate <- lapply(tx_LOHrate$total_tLOH, 
                                   function(x) sapply(chrom_lengths_BY_df$chrom_length,
                                                      function(y) x * (y/g_length))) %>% 
  as.data.frame(col.names = Tx_ID_levels) %>% 
  mutate(CHROM = chrom_lengths_BY_df$CHROM) %>% 
  pivot_longer(cols = 1:3, names_to = "Tx_ID", values_to = "expected_tLOH") %>%
  mutate(Tx_ID = factor(Tx_ID, Tx_ID_levels)) %>%
  arrange(Tx_ID, CHROM) %>% 
  select(expected_tLOH) %>%
  cbind(expected_CHROM_type_rate, .)

LOHcounts_CHROM_type_rate <- LOHcounts_CHROM_type_mean %>% 
  select(Tx_ID, CHROM, type, bp_rate) %>% 
  pivot_wider(names_from = "type", values_from = "bp_rate") %>%
  rename(bp_rate_iLOH = Interstitial, bp_rate_tLOH = Terminal) %>%
  arrange(Tx_ID, CHROM)

expected_CHROM_type_rate <- lapply(mean_LOHrate$sum_LOH,
                                function(x) sapply(chrom_lengths_BY_df$chrom_length,
                                                   function(y) x * (y/g_length) / y / n_gens /
                                                     tx_LOHrate$f_iLOH[tx_LOHrate$total_LOH == x] /
                                                     mean_LOHrate$n_clones[mean_LOHrate$sum_LOH == x])) %>%
  as.data.frame(col.names = Tx_ID_levels) %>%
  mutate(CHROM = chrom_lengths_BY_df$CHROM) %>%
  pivot_longer(cols = 1:3, names_to = "Tx_ID", values_to = "e_iLOH_bp_rate") %>%
  arrange(Tx_ID, CHROM) %>%
  select(e_iLOH_bp_rate) %>%
  cbind(expected_CHROM_type_rate, .)

expected_CHROM_type_rate$chrom_length <- rep(chrom_lengths_BY, 3)

expected_CHROM_type_long <- expected_CHROM_type_rate %>%
  ungroup() %>%
  rename(total_Interstitial = total_iLOH,
         total_Terminal = total_tLOH,
         expected_Interstitial = expected_iLOH,
         expected_Terminal = expected_tLOH) %>%
  pivot_longer(cols = !c(Tx_ID, CHROM, chrom_length), 
               names_to = "type", values_to = "n_LOH") %>%
  mutate(data_type = colsplit(type, "_", c("a", "b"))$a,
         type = colsplit(type, "_", c("a", "b"))$b)
expected_CHROM_type_long$type <- factor(expected_CHROM_type_long$type)
expected_CHROM_type_long$data_type <- factor(expected_CHROM_type_long$data_type)

expected_CHROM_type_long$clone_rate <- operate_by_factor_match(n_clones_LOH_xTx, 
                                                            expected_CHROM_type_long[c("Tx_ID", "n_LOH")],
                                                            function(x, y) y / x)
expected_CHROM_type_long$bp_rate <- expected_CHROM_type_long$clone_rate / 
  expected_CHROM_type_long$chrom_length /
  n_gens

chisq_test(expected_CHROM_type_long %>% 
             filter(Tx_ID == "W", type == "Interstitial") %>% 
             select("data_type", "n_LOH"))

rate_lm <- lm(bp_rate ~ chrom_length + Tx_ID + data_type, 
              data = expected_CHROM_type_long %>% filter(type == "Interstitial",
                                                         ))

summary(rate_lm)
# expected_TxCHROM_rate$per_gen <- expected_TxCHROM_rate$per_clone/n_gens

expected_CHROM_type_long %>% 
  filter(data_type == "total") %>%
  ggplot(aes(x = chrom_length, y = bp_rate, color = Tx_ID)) +
  geom_point(aes()) +
  geom_line(data = expected_CHROM_type_long %>% 
              filter(data_type == "expected"),
            aes(x = chrom_length, y = bp_rate, color = Tx_ID)) +
  scale_color_manual(values = txPal) +
  facet_grid(type~.)


###############################################################################
# Position specific effects of iLOH and tLOH ##########

POSi_data_in <- all_LOHbounds_merge_NS %>% get_LOH_coordinates()

POSi_het <- all_GT_bounds_merge %>% 
  filter(GT == "0/1")

POSi_data_wHet <- merge(POSi_data_in, POSi_het, all = T) %>% arrange(ID, est_start)

LOHcounts_Tx_type_mean <- POSi_data_in %>% group_by(Tx_ID, isTerm) %>% summarize(total_LOH = n())
LOHcounts_Tx_type_mean$per_clone <- operate_by_factor_match(n_clones_LOH_xTx, 
                                                            LOHcounts_Tx_type_mean[, c("Tx_ID", "total_LOH")],
                                                            .fun = function(x, y) y/x)

LOHcounts_ID_CHROM_type <- POSi_data_in %>% 
  filter(!isTerm) %>%
  group_by(ID, CHROM, .drop = F) %>% 
  count(name = "n_LOH") %>% mutate(type = "Interstitial") %>% as.data.frame()
LOHcounts_ID_CHROM_t <- POSi_data_in %>% 
  filter(isTerm) %>%
  group_by(ID, CHROM, .drop = F) %>%
  count(name = "n_LOH") %>% mutate(type = "Terminal") %>% as.data.frame()
LOHcounts_ID_CHROM_type <- rbind(LOHcounts_ID_CHROM_type, LOHcounts_ID_CHROM_t)

LOHcounts_ID_CHROM_type <- CategoriesFromID(LOHcounts_ID_CHROM_type)

LOHcounts_ID_CHROM_type$type <- factor(LOHcounts_ID_CHROM_type$type)
all_CHROM_type_LOHrate_perm <- data.frame(NULL)
for(ty in levels(LOHcounts_ID_CHROM_type$type)) {
  # ty = "Interstitial"
  for(tx in 1:nrow(Tx_ID_combo)){
    # tx = 1
    all_Chr_perm <- data.frame(Combo = with(Tx_ID_combo[tx, ], paste(Tx_1, Tx_2, sep = "-")),
                               type = ty,
                               CHROM = chrom_bound_BY$CHROM, 
                               obsvStat = 0, critVal = 0, pVal = 0, rejectNull = 0)
    for(c in seq_along(chrom_bound_BY$CHROM)) {
      # c = 14
      chr = chrom_bound_BY$CHROM[c]
      # Chr_LOHcounts_ID <- LOHcounts_ID_CHROM %>% filter(CHROM == chr)
      Chr_LOHcounts_ID <- LOHcounts_ID_CHROM_type %>% filter(CHROM == chr, type == ty)
      Chr_perm <- perm_test(Chr_LOHcounts_ID, 
                            cat_var = "Tx_ID", 
                            cat_names = with(Tx_ID_combo[tx, ], c(Tx_1, Tx_2)), 
                            response_var = "n_LOH", 
                            alpha = 0.05,
                            n_perms = 10000, alt_hyp = "two-tailed",
                            rtrn = "all", include_matrix = F)
      all_Chr_perm[c, 4:7] <- unlist(Chr_perm[1:4])
    }
    all_CHROM_type_LOHrate_perm <- rbind(all_CHROM_type_LOHrate_perm, all_Chr_perm)
  }
}

# all_CHROM_LOHrate_perm %>% arrange(desc(obsvStat)) %>% head()
# all_CHROM_LOHrate_perm %>% arrange(pVal) %>% head()

all_CHROM_type_LOHrate_slim <- all_CHROM_type_LOHrate_perm %>% filter(Combo != "C-D")
all_CHROM_type_LOHrate_slim <- BHcorrection(all_CHROM_type_LOHrate_slim)

###############################################################################
# Main Fig: Chromosome rates by Tx !!!!########################################

LOHcounts_CHROM_type_mean$Tx_ID <- Recode_Tx_ID(LOHcounts_CHROM_type_mean$Tx_ID, "Tx_ID")
LOHcounts_CHROM_type_mean$n_CHROM <- as.numeric(LOHcounts_CHROM_type_mean$CHROM)
LOHcounts_CHROM_type_mean$Tx_n <- as.numeric(LOHcounts_CHROM_type_mean$Tx_ID)

point_diff <- 0.16

LOHrate_type_CHROM_line <- LOHcounts_CHROM_type_mean %>% 
  ggplot() + 
  annotate(geom = "rect", xmin = LOHcounts_CHROM_type_mean$n_CHROM[c(F, T)] - 0.5, 
           xmax = LOHcounts_CHROM_type_mean$n_CHROM[c(F, T)] + 0.5, ymin = 0, ymax = Inf, 
           fill = "grey90", alpha = 0.1) +
  # geom_errorbar(aes(x = Tx_name, ymin = mean_rate, ymax = mean_rate, color = Tx_name), 
  #               width = 0.6, position = position_dodge(width = 0.6), size = 0.75) +
  geom_line(aes(x = n_CHROM + (Tx_n - 2) * point_diff, 
                y = bp_rate, color = Tx_ID, group = Tx_ID), 
            alpha = 0.7, size = 1) +
  geom_errorbar(aes(x = n_CHROM + (Tx_n - 2) * point_diff, 
                    ymin = rate_CI_lo, 
                    ymax = rate_CI_up, color = Tx_ID),
                width = 0, alpha = 0.7, size = 0.5) +
  geom_point(aes(x = n_CHROM + (Tx_n - 2) * point_diff, y = bp_rate, color = Tx_ID), 
             size = 4) +
  scale_color_manual(values = txPal, name = "Strain") +
  scale_y_continuous(labels = format_sci_10,
                     breaks = c(0, seq(1e-10, 5e-10, 1e-10))) +
  scale_x_continuous(labels = roman_chr, expand = c(0, 0),
                     breaks = 1:16, limits = c(0.5, 16.5)) +
  # ylim(c(0, NA)) +
  ylab("LOH event rate (/bp/gen)") + 
  xlab("Chromosome") +
  facet_wrap(type~., ncol = 1, strip.position = "top", scales = "free_y") +
  theme(legend.text = element_text(size = 16),
        legend.title = element_text(size = 20),
        strip.text = element_text(size = 24),
        axis.text.x = element_text(vjust = 2.5),
        axis.title.y = element_text(size = 24, vjust = 2.5),
        axis.title.x = element_text(size = 24),
        axis.text = element_text(size = 18),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        plot.margin = margin(20, 20, 20, 20),
        # panel.background = element_rect(color = "grey80"),
        # text = element_text(size = 24),
        legend.position = c(0.96, 0.94),
        legend.background = element_rect(fill = "white", color = "white"))

LOHrate_type_CHROM_line

ggsave(file.path(outIntDir, "LOHrate_CHROM_type_line_2022_04.png"), 
       plot = LOHrate_type_CHROM_line,
       device = "png",
       width = 16, height = 9, 
       units = "in",
       dpi = 600)

###############################################################################
# Empirical cumulative distribution of LOH rates ##############################
# along the length of each chromosome with permuted Anderson-Darling test
POSi_test_in <- POSi_data_in 

all_CHROM_type_LOHrate_perm <- data.frame(NULL)
for(ty in levels(LOHcounts_ID_CHROM_type$type)) {
  # ty = "Interstitial"
  for(tx in 1:nrow(Tx_ID_combo)){
    # tx = 1
    all_Chr_perm <- data.frame(Combo = with(Tx_ID_combo[tx, ], paste(Tx_1, Tx_2, sep = "-")),
                               type = ty,
                               CHROM = chrom_bound_BY$CHROM, 
                               obsvStat = 0, critVal = 0, pVal = 0, rejectNull = 0)
    for(c in seq_along(chrom_bound_BY$CHROM)) {
      # c = 14
      chr = chrom_bound_BY$CHROM[c]
      # Chr_LOHcounts_ID <- LOHcounts_ID_CHROM %>% filter(CHROM == chr)
      Chr_LOHcounts_ID <- POSi_test_in %>% filter(CHROM == chr, type == ty)
      
      Chr_perm <- perm_test(Chr_LOHcounts_ID, 
                            cat_var = "Tx_ID", 
                            cat_names = with(Tx_ID_combo[tx, ], c(Tx_1, Tx_2)), 
                            response_var = "n_LOH", 
                            alpha = 0.05,
                            n_perms = 10000, alt_hyp = "two-tailed",
                            rtrn = "all", include_matrix = F)
      all_Chr_perm[c, 4:7] <- unlist(Chr_perm[1:4])
    }
    all_CHROM_type_LOHrate_perm <- rbind(all_CHROM_type_LOHrate_perm, all_Chr_perm)
  }
}


all_CHROM_LOHprop_AD <- data.frame(NULL)
for(tx in 1:nrow(Tx_ID_combo)){
  # tx = 1
  Tx_1 <- Tx_ID_combo$Tx_1[tx]
  Tx_2 <- Tx_ID_combo$Tx_2[tx]
  Combo_Chr_AD <- data.frame(Combo = with(Tx_ID_combo[tx, ], paste(Tx_1, Tx_2, sep = "-")),
                             Tx_1 = Tx_1, Tx_2 = Tx_2,
                             CHROM = chrom_bound_BY$CHROM, 
                             ADstat = 0, pVal = 0)
  AD_LOHcounts <- POSi_test_in %>% filter(Tx_ID %in% c(Tx_1, Tx_2)) %>% 
    # filter(!isTerm) %>%
    select(Tx_ID, ID, CHROM, est_mid_POS)
  AD_LOHcounts$Tx_ID <- droplevels(AD_LOHcounts$Tx_ID)
  for(ch in seq_along(chrom_bound_BY$CHROM)) {
    # ch = 1
    chr = chrom_bound_BY$CHROM[ch]
    Chr_LOHcounts_ID <- AD_LOHcounts %>% filter(CHROM == chr)
    Chr_AD <- ad.test(est_mid_POS ~ Tx_ID, data = Chr_LOHcounts_ID, method = "exact")
    Combo_Chr_AD[ch, c(5, 6)] <- Chr_AD$ad[2, c(1, 3)]
  }
  all_CHROM_LOHprop_AD <- rbind(all_CHROM_LOHprop_AD, Combo_Chr_AD)
}

all_CHROM_LOHprop_AD_slim <- all_CHROM_LOHprop_AD %>% filter(Combo != "C-D")
all_CHROM_type_LOHprop_AD_slim <- rbind(all_CHROM_LOHprop_AD_slim,
                                        all_CHROM_LOHprop_AD %>% 
                                          filter(Combo != "C-D") %>% mutate(type = "Terminal"))
all_CHROM_type_LOHprop_AD_slim <- BHcorrection(all_CHROM_type_LOHprop_AD_slim)



LOHcounts_ID_CHROM_type$chrom_length <- 
  operate_by_factor_match(chrom_lengths_BY_df[, c("CHROM", "chrom_length")],
                          LOHcounts_ID_CHROM_type[, c("CHROM", "Tx_ID")],
                          function(x, y) x)


###############################################################################
# LOH rate vs chromsome size ####


lm(n_LOH ~ CHROM_len*type, 
   data = LOHcounts_ID_CHROM_type) %>% summary()

###############################################################################
# Breakpoint /bp/clone rate vs position #####
iLOH_BPrate_SW <- POSi_data_in %>% 
  mutate(bp = 1) %>%
  filter(!isTerm) %>%
  # filter(!cross_h_term) %>% 
  # mutate(POS = fract_dist) %>%
  mutate(POSi = est_mid) %>%
  SliderCalc(data_col = "bp", index_col = "POSi", 
             window_size = 50000, slide_interval = 5000,
             summary_stat = sum, factor_col = "Tx_ID", chrom_win = T)
iLOH_BPrate_SW$Tx_ID <- factor(iLOH_BPrate_SW$Tx_ID, levels = Tx_ID_levels)

iLOH_BPrate_SW$start_POS <- iLOH_BPrate_SW %>% ConvertPosIndicies(pos_col = "start", index_out = "POS")
iLOH_BPrate_SW$end_POS <- iLOH_BPrate_SW %>% ConvertPosIndicies(pos_col = "end", index_out = "POS")
iLOH_BPrate_SW$mid_POS <- round((iLOH_BPrate_SW$start_POS + iLOH_BPrate_SW$end_POS)/2)
iLOH_BPrate_SW$length <- iLOH_BPrate_SW$end - iLOH_BPrate_SW$start + 1
iLOH_BPrate_SW$per_clone <- operate_by_factor_match(n_clones_LOH_xTx, 
                                                   iLOH_BPrate_SW[, c("Tx_ID", "n_elements")],
                                                   .fun = function(x, y) y/x)

iLOH_BPrate_SW$type <- "Interstitial"

tLOH_BPrate_SW <- POSi_data_in %>% 
  mutate(bp = 1) %>%
  filter(isTerm) %>%
  # filter(!cross_h_term) %>% 
  # mutate(POS = fract_dist) %>%
  mutate(POSi = est_mid) %>%
  SliderCalc(data_col = "bp", index_col = "POSi", 
             window_size = 50000, slide_interval = 5000,
             summary_stat = sum, factor_col = "Tx_ID", chrom_win = T)
tLOH_BPrate_SW$Tx_ID <- factor(tLOH_BPrate_SW$Tx_ID, levels = Tx_ID_levels)

tLOH_BPrate_SW$start_POS <- tLOH_BPrate_SW %>% ConvertPosIndicies(pos_col = "start", index_out = "POS")
tLOH_BPrate_SW$end_POS <- tLOH_BPrate_SW %>% ConvertPosIndicies(pos_col = "end", index_out = "POS")
tLOH_BPrate_SW$mid_POS <- round((tLOH_BPrate_SW$start_POS + tLOH_BPrate_SW$end_POS)/2)
tLOH_BPrate_SW$length <- tLOH_BPrate_SW$end - tLOH_BPrate_SW$start + 1
tLOH_BPrate_SW$per_clone <- operate_by_factor_match(n_clones_LOH_xTx, 
                                                    tLOH_BPrate_SW[, c("Tx_ID", "n_elements")],
                                                    .fun = function(x, y) y/x)
tLOH_BPrate_SW$type <- "Terminal"


bLOH_BPrate_SW <- rbind(iLOH_BPrate_SW, tLOH_BPrate_SW)
bLOH_BPrate_SW$type <- factor(bLOH_BPrate_SW$type)
bLOH_BPrate_SW <- bLOH_BPrate_SW %>% 
  mutate(bp_rate = per_clone/length/n_gens)

bLOH_BPrate_SW <- bLOH_BPrate_SW %>% 
  mutate(mid_POSi = (start + end)/2,
         mid_POS_kb = mid_POS/1000, mid_POSi_kb = (start + end)/2000)

bLOH_BPrate_SW$rom_CHROM <- factor(chrom_IDs$rom_CHROM[as.numeric(bLOH_BPrate_SW$CHROM)], 
                                   levels = roman_chr)


max_rate <- max(bLOH_BPrate_SW$bp_rate) * 1.2
y_scale <- ceiling(log10(max_rate)) - 1
y_max <- 10^y_scale
y_incr <- 2.5*10^(y_scale - 1)

BPrate_iLOH_SW_plot <- bLOH_BPrate_SW %>% 
  filter(type == "Interstitial") %>%
  filter(bp_rate < 1E-8) %>%
  # filter(length > 10000) %>%
  ggplot(aes(x = start_POS/1000, y = bp_rate)) + 
  geom_vline(data = centrom_df, aes(xintercept = POS/1000), 
             linetype = 2,color = "blue4", size = 0.25) +
  # geom_ribbon(aes(ymin = CI_lo, ymax = CI_up, fill = Tx_ID, color = Tx_ID), 
  #             size = 0.05, alpha = 0.1) +
  geom_line(aes(color = Tx_ID), size = 1) +
  # geom_line(data = genome_scores_50kb, 
  #           aes(x = mid_POS/1000, y = -(score - min_score)/((max(score) - min_score)*3))) +
  # geom_line(aes(x = mid_POS_kb, y = log(rate, 10), color = Tx_ID)) +
  # geom_point(aes(color = Tx_ID), size = 0.5) +
  # xlim(0, max(chrom_lengths_BY)) +
  xlab("Window Position (kb)") + xlim(0, NA) +
  scale_y_continuous(labels = format_sci_10,
                     # labels = c(0, format(seq(y_incr, y_max, y_incr), scientific = T)),
                     breaks = c(0, seq(y_incr, y_max, y_incr)),
                     name = "Breakpoint rate /bp/gen") +
  scale_color_manual(values = txPal, name = "Strain") + 
  scale_fill_manual(values = txPal, name = "Strain") + 
  facet_wrap(~rom_CHROM, ncol = 4, scales = "free_x") +
  theme(legend.text = element_text(size = 16),
        legend.title = element_text(size = 18),
        strip.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 16),
        # text = element_text(size = 18), 
        plot.margin = margin(20, 20, 20, 20),
        panel.grid.minor = element_blank(),
        legend.position = c(0.96, 0.95),
        legend.background = element_rect(fill = "white", color = "white"))

BPrate_iLOH_SW_plot


BPrate_tLOH_SW_plot <- bLOH_BPrate_SW %>% 
  filter(type == "Terminal") %>%
  filter(bp_rate < 7.5E-9) %>%
  # filter(length > 10000) %>%
  ggplot(aes(x = start_POS/1000, y = bp_rate)) + 
  geom_vline(data = centrom_df, aes(xintercept = POS/1000), 
             linetype = 2,color = "blue4", size = 0.25) +
  # geom_ribbon(aes(ymin = CI_lo, ymax = CI_up, fill = Tx_ID, color = Tx_ID), 
  #             size = 0.05, alpha = 0.1) +
  geom_line(aes(color = Tx_ID), size = 1) +
  # geom_line(data = genome_scores_50kb, 
  #           aes(x = mid_POS/1000, y = -(score - min_score)/((max(score) - min_score)*3))) +
  # geom_line(aes(x = mid_POS_kb, y = log(rate, 10), color = Tx_ID)) +
  # geom_point(aes(color = Tx_ID), size = 0.5) +
  # xlim(0, max(chrom_lengths_BY)) +
  xlab("Window Position (kb)") + xlim(0, NA) +
  scale_y_continuous(labels = format_sci_10,
                     # labels = c(0, format(seq(y_incr, y_max, y_incr), scientific = T)),
                     breaks = c(0, seq(y_incr, y_max, y_incr)),
                     name = "Breakpoint rate /bp/gen") +
  scale_color_manual(values = txPal, name = "Strain") + 
  scale_fill_manual(values = txPal, name = "Strain") + 
  facet_wrap(~rom_CHROM, ncol = 4, scales = "free_x") +
  theme(strip.text = element_text(size = 20),
        axis.title.y = element_text(size = 24, vjust = 2.5),
        axis.title.x = element_text(size = 24, vjust = -1.7),
        axis.text = element_text(size = 18),
        # text = element_text(size = 18), 
        plot.margin = margin(20, 20, 20, 20),
        panel.grid.minor = element_blank(),
        legend.position = "none")


BPrate_tLOH_SW_plot

ggsave(file.path(outIntDir, "BPrate_SW_SW50x5_plot_2022_04.png"), 
       plot = BPrate_SW_plot,
       device = "png",
       width = 16, height = 9, 
       units = "in",
       dpi = 600)



figure <- plot_grid(BPrate_iLOH_SW_plot, BPrate_tLOH_SW_plot,
                    labels = c("A", "B"),
                    align = "v",
                    scale = 0.95,
                    label_size = 28,
                    hjust = 0,
                    ncol = 1, nrow = 2) +
  theme(plot.background = element_rect(fill = "white", colour = "white"))

figure

ggsave(file.path(outIntDir, "SNM_Indel_combo_2022_05.png"), 
       plot = figure,
       device = "png",
       width = 11, height = 11, 
       units = "in",
       dpi = 600)



POSi_data_wHet$GT <- factor(POSi_data_wHet$GT, levels = GT_levels)