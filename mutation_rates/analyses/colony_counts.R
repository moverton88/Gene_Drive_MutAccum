# To test the degree of mutational burden in the Drive treatment, 
# XC, SG, and AM revived one ancestor and all decendant end-point clones 
# from the freezer and measured colony sizes after growth under the MA conditions
# These data were logged in a spreadsheet here:
# https://docs.google.com/spreadsheets/d/1T7h2aQnJkv2Ok9ng8TLzKfqB3Z7zCPwaE3K1gPicy_M/edit?usp=sharing
# though I reformatted the data for import into R here:
# https://docs.google.com/spreadsheets/d/1T7h2aQnJkv2Ok9ng8TLzKfqB3Z7zCPwaE3K1gPicy_M/edit?usp=sharing

col_classes <- c("ffccfiiidld")
colony_counts <- read_delim(colony_fn, delim = "\t", col_types = col_classes) %>% as.data.frame()
colnames(colony_counts) <- c("Tx_name", "clone", "clone_ID", "sample_ID", 
                             "person", "rep", "e_3", "e_4", "ratio", 
                             "valid", "count_e_3")

colony_counts$clone_ID <- factor(str_pad(colony_counts$clone_ID, 3, pad = "0"))

colony_counts_valid <- colony_counts %>% filter(valid)

colony_counts_valid <- colony_counts_valid %>% mutate(count_comb = round((e_3 + e_4)/1.1*1000))

colony_counts_valid %>% mutate(diff = (count_e_3 - count_comb)/count_e_3) %>% arrange(diff) %>% head()

count_stats <- colony_counts_valid %>% group_by(Tx_name, clone_ID) %>% 
  summarize(mean_count = mean(count_comb, na.rm = T), 
            stdev = sd(count_comb, na.rm = T), 
            sd_over_pois = sd(count_comb, na.rm = T)/sqrt(mean(count_comb, na.rm = T)),
            range_over_mean = (max(count_comb, na.rm = T) - min(count_comb, na.rm = T))/mean(count_comb, na.rm = T),
            n_reps = n())

colony_counts_valid %>% 
  ggplot() + 
  geom_point(aes(x = clone_ID, y = count_comb, color = Tx_name), size = 4, alpha = 0.6) + 
  scale_color_manual(values = txPal) +
  ylab("estimated cell count") + scale_y_continuous(labels = c(0, format(seq(5e5, 1.5e6, 5e5), scientific = T)),
                                               breaks = c(0, seq(5e5, 1.5e6, 5e5)))

count_stats %>% 
  filter(n_reps > 1) %>% 
  ggplot() + 
  geom_point(aes(x = clone_ID, y = mean_count, color = Tx_name), size = 2) + 
  geom_errorbar(aes(x = clone_ID, y = mean_count, 
                    ymin = mean_count - stdev, ymax = mean_count + stdev, 
                    color = Tx_name), width = 0.5, size = 0.7) +
  scale_color_manual(values = txPal, name = "Strain") +
  ylab("mean count (+/- sd)")

nb_model <- data.frame()

lm_var <- lm(I(stdev^2) ~ 0 + I(mean_count), data = count_stats)
summary(lm_var)
lm_b_var <- lm_var$coefficients[1]
lm_b_var <- 0
lm_m_var <- lm_var$coefficients[1]
lm_m_var <- 0
lm_q_var <- lm_var$coefficients[1]
lm_q_var <- 0
x_max <- max(count_stats$mean_count)
x_seg <- round(seq.int(0, x_max, length.out = 30))
x_seg_end <- x_seg[2:length(x_seg)]
x_seg <- x_seg[1:(length(x_seg) - 1)]
lm_segment <- data.frame(x = x_seg, 
                         y = lm_b_var + lm_m_var*x_seg + lm_q_var*x_seg^2, 
                         x_end = x_seg_end, 
                         y_end = lm_b_var + lm_m_var*x_seg_end + lm_q_var*x_seg_end^2)

count_stats %>% 
  filter(n_reps > 2) %>% 
  # filter(stdev < 200000) %>%
  ggplot() + 
  geom_abline(color = "grey30") + 
  geom_segment(data = lm_segment, aes(x = x, y = y, xend = x_end, yend = y_end)) +
  geom_point(aes(x = mean_count, y = stdev^2, color = Tx_name), size = 3) + 
  scale_color_manual(values = txPal) 




anc_counts <- colony_counts_valid %>% filter(clone == "Anc")
anc_count_means <- anc_counts %>% group_by(Tx_name) %>% 
  summarize(mean_count = mean(count_e_3), mean_count_comb = mean(count_comb))
end_counts <- colony_counts_valid %>% filter(clone == "End")
end_count_means <- end_counts %>% group_by(Tx_name, clone_ID) %>% 
  summarize(mean_count = mean(count_e_3), mean_count_comb = mean(count_comb), mean_nrml_cc = mean(nrml_count_comb))

end_counts$nrml_count_comb <- 0
for(tx in Tx_name_levels) {
  # tx <- "WT"
  end_counts$nrml_count_comb[end_counts$Tx_name == tx] <- end_counts$count_comb[end_counts$Tx_name == tx] / 
    anc_count_means$mean_count_comb[anc_count_means$Tx_name == tx]
}

end_count_means %>% ggplot(aes(x = Tx_name, y = mean_nrml_cc)) + geom_point() + stat_summary(geom = "point", size = 5, fun = mean)

colony_aov <- aov(nrml_count_comb ~ Tx_name, data = end_counts)
summary(colony_aov)

colony_wrs <- wilcox.test(mean_nrml_cc ~ Tx_name, data = end_count_means %>% filter(Tx_name != "Drive"), distribution = "exact")
colony_wrs
