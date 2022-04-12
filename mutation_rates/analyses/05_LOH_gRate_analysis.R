## Estimate and compare the number of LOH events
# Uses the function "EstDataBounds" to estimate the LOH regions
# Then the "MarkLOHcomplexes" to mark LOH regions within recombination
# tracts
# Then the "CountLOHevents" function to count simple and complex events
# These three functions are located in the LOHanalysisFunctions.R file
# The "perm_test" from the file permutationFunction.R is used to test 
# for differences in means among treatments and to test whether the 
# LOH number distributions follow a poisson distribution

# Load in counts table
# LOHcounts_in <- all_LOHcounts_merge
LOHcounts_in <- all_LOHcounts_merge_NS
# LOHcounts_in <- all_LOHcounts_merge_HE
# LOHcounts_in <- all_LOHcounts_merge_EC

LOHcounts_in$Tx_ID <- Recode_Tx_ID(LOHcounts_in$Tx)
LOHcounts_in <- LOHcounts_in %>% arrange(Tx_ID, ID)

# LOHcounts_in %>% arrange(desc(n_LOH)) %>% head()

# Permutation test of equal median n_LOH for ea vs ea Tx --------
LOHcounts_in %>% group_by(Tx_ID) %>% summarize(median_nLOH = median(n_LOH))
LOHcounts_in %>% # filter(!Line %in% c("H_F", "H_H")) %>% 
  group_by(Tx_ID) %>% summarize(mean_nLOH = mean(n_LOH))


LOH_final_counts_stats <- LOHcounts_in %>% 
  group_by(Tx_ID) %>% 
  summarize(n_clones = n(),
            tot_LOH = sum(n_LOH),
            mean_n_LOH = mean(n_LOH), 
            # se_n = se(n_LOH), 
            meanCI95lo = (mean(n_LOH) - se(n_LOH)*1.96), 
            meanCI95hi = (mean(n_LOH) + se(n_LOH)*1.96),
            genrate = mean(n_LOH)/n_gens, 
            genrateCI95lo = (mean(n_LOH) - se(n_LOH)*1.96)/n_gens, 
            genrateCI95hi = (mean(n_LOH) + se(n_LOH)*1.96)/n_gens)

W_C_mean_perm <- LOHcounts_in %>% 
  # filter(!Line %in% c("H_F", "H_H")) %>%
  perm_test(df=., cat_var ="Tx_ID", cat_names = c("C", "W"), response_var = "n_LOH", 
            n_perms = 10000, alpha = 0.01, alt_hyp = "two-tailed", test_stat = mean)

W_D_mean_perm <- LOHcounts_in %>% 
  perm_test(df=., cat_var ="Tx_ID", cat_names = c("D", "W"), response_var = "n_LOH", 
            n_perms = 10000, alpha = 0.01, alt_hyp = "two-tailed", test_stat = mean)

C_D_mean_perm <- LOHcounts_in %>% 
  perm_test(df=., cat_var ="Tx_ID", cat_names = c("D", "C"), response_var = "n_LOH", 
            n_perms = 10000, alpha = 0.01, alt_hyp = "two-tailed", test_stat = mean)

# Permutation test of equal variance across Tx -------------
W_C_IQR_perm <- LOHcounts_in %>% 
  # filter(!Line %in% c("H_F", "H_H")) %>%
  perm_test(df=., cat_var ="Tx_ID", cat_names = c("C", "W"), response_var = "n_LOH", 
            n_perms = 10000, alpha = 0.01, alt_hyp = "two-tailed", test_stat = IQR)

W_D_IQR_perm <- LOHcounts_in %>% 
  perm_test(df=., cat_var ="Tx_ID", cat_names = c("D", "W"), response_var = "n_LOH", 
            n_perms = 10000, alpha = 0.01, alt_hyp = "two-tailed", test_stat = IQR)

C_D_IQR_perm <- LOHcounts_in %>% 
  perm_test(df=., cat_var ="Tx_ID", cat_names = c("D", "C"), response_var = "n_LOH", 
            n_perms = 10000, alpha = 0.01, alt_hyp = "two-tailed", test_stat = IQR)

mean_p <- c(round(W_D_mean_perm$p_value, 3),
            round(W_C_mean_perm$p_value, 3),
            round(C_D_mean_perm$p_value, 3))
median_is_sig <- ifelse(mean_p <= 0.01, "*", "")

IQR_p <- c(round(W_D_IQR_perm$p_value, 3),
           round(W_C_IQR_perm$p_value, 3),
           round(C_D_IQR_perm$p_value, 3))
IQR_is_sig <- ifelse(IQR_p <= 0.05, "*", "")

W_D_test_label <- paste0("mean p = ", mean_p[1], median_is_sig[1],
                         "\nIQR p = ", IQR_p[1], IQR_is_sig[1])
W_C_test_label <- paste0("mean p = ", mean_p[2], median_is_sig[2],
                         "\nIQR p = ", IQR_p[2], IQR_is_sig[2])
C_D_test_label <- paste0("mean p = ", mean_p[3], median_is_sig[3],
                         "\nIQR p = ", IQR_p[3], IQR_is_sig[3])


LOHcounts_in$dot_color <- txPal[as.numeric(LOHcounts_in$Tx_name)]
LOHcounts_in$fill_color <- txPal[as.numeric(LOHcounts_in$Tx_name)]
LOHcounts_in$fill_color[LOHcounts_in$ID %in% DA_clones] <- NA
LOHcounts_in <- LOHcounts_in %>% arrange(Tx_ID, n_LOH)
LOH_final_counts_stats$fill_color <- rev(txPal[as.numeric(LOH_final_counts_stats$Tx_ID)])

lab_y <- max(LOHcounts_in$n_LOH)
nudge <- c(2, 6, 10)
max_break <- ceiling(lab_y/5)*5
sr <- 1.1
x_diff <- (sr - 1)/50

LOHcount_dotplot_mean <- LOHcounts_in %>%
  # filter(Tx_name == "Cas9") %>%
  ggplot(aes(x = x_diff)) + 
  geom_dotplot(aes(y = n_LOH, group = Tx_ID), 
               fill = LOHcounts_in$fill_color, 
               color = LOHcounts_in$dot_color,
               stroke = 2,
               binwidth = 1, binaxis = "y", 
               stackdir = "center", binpositions = "all",
               dotsize = 0.4, stackratio = sr) +
  stat_summary(aes(y = n_LOH), color = "grey30", fun = "mean", fun.min = "mean", fun.max= "mean",
               size = 0.2, width = 0.9, geom = "crossbar") +
  scale_fill_manual(values = txPal) +
  # scale_color_manual(values = txPal) +
  scale_y_continuous(breaks = seq(0, max_break, 5),
                     name = "Number of LOH events") +
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
        text = element_text(size = 30), 
        strip.text = element_text(size = 30),
        # axis.text = element_text(size = 12),
        legend.position = "none") +
  # facet_grid(.~Line %in% c("H_F", "H_H"))
  facet_grid(.~Tx_ID)

LOHcount_dotplot_mean

ggsave(file.path(outIntDir, "LOHcount_dotplot_2022_03_DE.png"),
       plot = LOHcount_dotplot_mean,
       device = "png",
       width = 16, height = 9,
       units = "in",
       dpi = 600)

figure <- plot_grid(LOHcount_dotplot_mean, SNMcount_dotplot_mean,
                    labels = c("A", "B"),
                    align = "h",
                    scale = 1,
                    label_size = 20,
                    hjust = 0,
                    ncol = 1, nrow = 2)

figure

ggsave(file.path(outIntDir, "MutCount_combo_2022_03.png"), 
       plot = figure,
       device = "png",
       width = 11, height = 8.5, 
       units = "in",
       dpi = 600)

# Show which Drive clones were tested for drive activity ######################
DA_clones <- c("F_A09", "F_C02", "F_D01", "F_F03", "F_F07", "F_G10")
DA_counts_in <- LOHcounts_in %>% filter(Tx_name == "Drive") %>% arrange(n_LOH) %>% mutate(DA = as.numeric(ID %in% DA_clones) * 0.7 + 0.3)
DA_dotplot_mean <- ggplot(data = DA_counts_in, aes(x = 0, y = n_LOH)) +
  geom_dotplot(alpha = DA_counts_in$DA, fill = txPal[3], 
               color = "white", binwidth=1, binaxis = "y", 
               stackdir = "center", binpositions = "all",
               dotsize = 0.5, stackratio = 1) + 
  # geom_bracket(xmin = c("WT", "WT", "Cas9"), xmax = c("Cas9", "Drive", "Drive"),
  #              y.position = c(lab_y + nudge[1], lab_y + nudge[2], lab_y + nudge[3]),
  #              label = c(N_H_test_label, N_F_test_label, H_F_test_label),
  #              tip.length = 0.025, label.size = 4) +
  # scale_y_continuous(expand = c(0,0), limit = c(0, lab_y + nudge[3] + 5), 
  #                    breaks = seq(0, max_break, 5),
  #                    name = "Number of LOH events") +
  # scale_alpha_discrete(range = c(0.3, 1)) +
  # scale_fill_manual(values = c("grey80", "red3")) +
  scale_y_continuous(breaks = seq(0, max_break, 5),
                     name = "Number of LOH events") +
  scale_x_continuous(name = "", breaks = 0, labels = "Drive") +
  theme(axis.ticks.y = element_blank(), 
        # axis.ticks.x = element_line(size = 0.25, color = "grey80"),
        # axis.ticks.length = unit(4, "mm"),
        plot.margin = unit(c(t = 5, r = 5, b = 5, l = 5), "mm"),
        # panel.grid.minor.x=element_blank(), 
        # panel.grid.major.x=element_line(),
        # panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        # axis.title.y = element_text(hjust = 0.17, vjust = 1.5),
        axis.title.x = element_text(vjust = -1.5),
        # axis.text.x = element_text(vjust = 1, margin = unit(c(t = 3, r = 0, b = 0, l = 0), "mm")),
        text = element_text(size = 18), 
        # axis.text = element_text(size = 12),
        legend.position = "none")

DA_dotplot_mean

ggsave(file.path(outIntDir, "DA_dotplot_mean_2022_02.png"), 
       plot = DA_dotplot_mean,
       device = "png",
       width = 8, height = 11.5, 
       units = "in",
       dpi = 600)

N_H_test_label <- paste0("p = ", median_p[1], median_is_sig[1])
N_F_test_label <- paste0("p = ", median_p[2], median_is_sig[2])
H_F_test_label <- paste0("p = ", median_p[3], median_is_sig[3])

nudge <- c(2, 4, 6)

LOHcount_dotplot_med <- LOHcounts_in %>%
  ggplot(aes(x = Tx_name)) + 
  annotate(geom = "rect", xmin = 0.9, xmax = 3.1, 
           ymin = lab_y + nudge[1] - 1, ymax = lab_y + nudge[3] + 5, 
           color = "white", fill = "white") +
  stat_summary(aes(y=n_LOH), color = "grey20", fun = "median", fun.min = "median", fun.max= "median",
               size= 0.5, width = 0.9, geom = "crossbar") +
  geom_dotplot(aes(y=n_LOH, fill = Tx_name), color = "white", binwidth=1, binaxis = "y", 
               stackdir = "center", binpositions = "all",
               dotsize = 0.5, stackratio = 1) + 
  geom_bracket(xmin = c("WT", "WT", "Cas9"), xmax = c("Cas9", "Drive", "Drive"),
               y.position = c(lab_y + nudge[1], lab_y + nudge[2], lab_y + nudge[3]),
               label = c(N_H_test_label, N_F_test_label, H_F_test_label),
               tip.length = 0.025, label.size = 5) +
  scale_fill_manual(values = txPal) +
  scale_color_manual(values = txPal) +
  scale_y_continuous(expand = c(0,0), limit = c(0, lab_y + nudge[3] + 5), breaks = seq(0, max_break, 5),
                     name = "Number of LOH events") +
  scale_x_discrete(name = "Strain") +
  theme(axis.ticks.y= element_blank(), 
        axis.ticks.x = element_line(size = 0.25, color = "grey80"),
        axis.ticks.length = unit(4, "mm"),
        plot.margin = unit(c(t = 5, r = 5, b = 5, l = 5), "mm"),
        panel.grid.minor.x=element_blank(), 
        panel.grid.major.x=element_line(),
        # panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        # axis.title.y = element_text(hjust = 0.17),
        axis.title.x = element_text(vjust = -1.5),
        axis.text.x = element_text(vjust = 1, margin = unit(c(t = 3, r = 0, b = 0, l = 0), "mm")),
        axis.title.y = element_text(hjust = 0.13, vjust = 1.7),
        text = element_text(size = 18), 
        axis.text = element_text(size = 16),
        legend.position = "none")

LOHcount_dotplot_med

ggsave(file.path(outIntDir, "LOHcount_dotplot_EC_2022_02.png"), 
       plot = LOHcount_dotplot_med,
       device = "png",
       width = 11, height = 8.5, 
       units = "in",
       dpi = 600)

# Proportion of genome converted ##############################################

prop_conv <- all_LOHbounds_merge %>% 
  filter(!is_error) %>% group_by(ID) %>% 
  summarize(p_conv = sum(est_length)/(mean_cover * 2))
prop_conv <- CategoriesFromID(prop_conv)
prop_conv_means <- prop_conv %>% group_by(Tx_name) %>% summarize(p = mean(p_conv))

prop_test_WTvDrive <- prop_conv %>% 
  perm_test(cat_var = "Tx_name", cat_names = c("WT", "Drive"), 
            response_var = "p_conv", test_stat = mean)

prop_test_WTvCas9 <- prop_conv %>% 
  perm_test(cat_var = "Tx_name", cat_names = c("WT", "Cas9"), 
            response_var = "p_conv", test_stat = mean)

prop_conv %>% ggplot() + 
  geom_histogram(aes(x = p_conv, y = ..density.., fill = Tx_name), center = 0, binwidth = 0.025) + 
  geom_vline(data = prop_conv_means, aes(xintercept = p)) +
  scale_fill_manual(values = txPal) +
  facet_grid(Tx_name~.)

prop_conv_aov <- aov(p_conv ~ Tx_name, data = prop_conv)
prop_conv_t_WTvDrive <- t.test(p_conv ~ Tx_name, data = prop_conv %>% filter(Tx_name != "Cas9"))
prop_conv_t_WTvCas9 <- t.test(p_conv ~ Tx_name, data = prop_conv %>% filter(Tx_name != "Drive"))
prop_conv_w_WTvDrive <- wilcox.test(p_conv ~ Tx_name, data = prop_conv %>% filter(Tx_name != "Cas9"))
prop_conv_w_WTvCas9 <- wilcox.test(p_conv ~ Tx_name, data = prop_conv %>% filter(Tx_name != "Drive"))

prop_conv %>% group_by(Tx_name) %>% summarize(med_p_conv = median(p_conv))

# Bootstrap confidence intervals ------
# Generation data is from https://drive.google.com/drive/folders/1vzY1GIyK2VaxDxHDtOYjV8vBSO139hEk?usp=sharing
mean.fun <- function(d, idx) {
  mean((d[idx]), trim = 0, na.rm = T)
}

sum.fun <- function(d, idx) {
  sum((d[idx]), trim = 0, na.rm = T)
}

LOH_CIs <- data.frame(NULL)
for(tx in Tx_name_levels) {
  tx_n_LOH <- LOHcounts_in %>% filter(Tx_name == tx) %>% pull(n_LOH)
  tx_boot <- boot(tx_n_LOH, statistic = mean.fun, R = 10000)
  tx_boot_CI <- boot.ci(boot.out = tx_boot,
                        type = "basic")$basic[4:5]
  tx_CI <- data.frame(Tx_name = tx, low95CI = tx_boot_CI[1], up95CI = tx_boot_CI[2])
  LOH_CIs <- rbind(LOH_CIs, tx_CI)
}

LOH_CIs$low95CI_rate <- LOH_CIs$low95CI / n_gens
LOH_CIs$up95CI_rate <- LOH_CIs$up95CI / n_gens

mean_LOHrate <- LOHcounts_in %>% 
  ungroup() %>% group_by(Tx_name) %>% 
  summarise(mean_LOH = mean(n_LOH), sd_LOH = sd(n_LOH), n_clones = n())

mean_LOHrate$se_LOH <- mean_LOHrate$sd_LOH/sqrt(mean_LOHrate$n_clones)
mean_LOHrate$mean_rate <- mean_LOHrate$mean_LOH / n_gens
mean_LOHrate$se_rate <- mean_LOHrate$se_LOH / n_gens

tx_LOHrate <- LOHcounts_in %>% 
  ungroup() %>% group_by(Tx_name) %>% 
  summarise(total_LOH = sum(n_LOH))
colnames(tx_LOHrate)[1] <- "Strain"

tx_LOHrate$n_clones <- LOHcounts_in %>% ungroup() %>% count(Tx_name) %>% pull(n)
tx_LOHrate$rate <- signif(tx_LOHrate$total_LOH / tx_LOHrate$n_clones / n_gens, 3)
tx_LOHrate$low95CI_rate <- LOH_CIs$low95CI_rate
tx_LOHrate$up95CI_rate <- LOH_CIs$up95CI_rate
tx_LOHrate <- as.data.frame(tx_LOHrate)

LOH_rate_table_out <- htmlTable(tx_LOHrate)

save_kable(LOH_rate_table_out, file = paste0(outIntDir, "LOH_rate_table_2022_02.png"))

# Output table as .png
png(paste0(outIntDir, "LOH_rate_table_2022_02.png"), 
    height = 50 * nrow(tx_LOHrate), width = 100 * ncol(tx_LOHrate))
grid.table(d = tx_LOHrate, theme = theme_minimal(), rows = NULL)
dev.off()

# Allele bias
allele_LOH_counts <- LOHcounts_in %>% summarize(ID = ID, Tx_name = Tx_name, n_BY_LOH = n_BYsmpl + n_BYcmplx, n_RM_LOH = n_RMsmpl + n_RMcmplx)
sum(allele_LOH_counts$n_BY_LOH)
sum(allele_LOH_counts$n_RM_LOH)

allele_LOH_counts %>% 
  group_by(Tx_name) %>% 
  summarise(n_BY = sum(n_BY_LOH), n_RM = sum(n_RM_LOH), bias = sum(n_BY_LOH)/(sum(n_BY_LOH)+sum(n_RM_LOH)))

###############################################################################
# Test whether LOH accumulation is from Poisson ----------

# Poisson test for pooled data
nRuns_freq <- LOHcounts_in %>% count(n_LOH)
mean_nRuns <- mean(LOHcounts_in$n_LOH)
pois_dist <- dpois(nRuns_freq$n_LOH, mean_nRuns, log = F)
resid <- 1-sum(pois_dist)
resid_prop <- (pois_dist+resid/length(pois_dist))*resid
pois_dist_crct <- pois_dist + resid_prop

chi_pois_test <- chisq.test(nRuns_freq$n, p=pois_dist_crct)

ggplot() + geom_col(aes(x=nRuns_freq$n_LOH, y=nRuns_freq$n/sum(nRuns_freq$n))) +
  geom_line(aes(x=nRuns_freq$n_LOH, y=pois_dist_crct), color="red") +
  geom_label(aes(x = 12, y = 0.18, label = paste0("X^2 p = ", signif(chi_pois_test$p.value, 3))), label.size = 0)

#  Poisson test for Tx data
nRuns_freq_Tx <- LOHcounts_in %>% count(Tx_name, n_LOH)
colnames(nRuns_freq_Tx) <- c("Tx_name", "value", "cnt")
nRuns_freq_Tx$Tx_name <- factor(nRuns_freq_Tx$Tx_name, levels = c("WT", "Cas9", "Drive"))
nRuns_freq_Tx <- nRuns_freq_Tx %>% group_by(Tx_name) %>% mutate(prop = cnt/sum(cnt))

# Generate Poisson dist for plotting 
pois_dist_tx_crct <- list()
for (t in 1:3) {
  # t = 1
  tx <- levels(nRuns_freq_Tx$Tx_name)[t]
  idx <- nRuns_freq_Tx$Tx_name == tx
  idx2 <- LOHcounts_in$Tx_name == tx
  pois_dist_tx <- dpois(0:max(nRuns_freq_Tx$value), mean(LOHcounts_in$n_LOH[idx2]), log = F)
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
pois_df$Tx_ID <- substr(pois_df$Tx_name, 1, 1)
pois_df$Tx_ID <- factor(pois_df$Tx_ID, levels = c("D", "C", "W"))

# Estimate p-values for null hyp of equal variances by Tx ------
LOHcount_pois <- LOHcounts_in %>% select(ID, n_LOH, Tx_ID) %>% group_by(Tx_ID) %>% 
  mutate(pois = rpois(n = n(), lambda = mean(n_LOH)))

perm_test(LOHcount_pois, cat_var = "Tx_ID", cat_names = c("W", "D"), response_var = 

nTrials <- 10000
var_list <- list()
var_sample <- list()
p_vals <- data.frame(NULL)
for (t in 1:3) {
  tx <- levels(LOHcounts_in$Tx_name)[t]
  idx2 <- LOHcounts_in$Tx_name == tx
  var_sample[t] <- var(LOHcounts_in$n_LOH[idx2])
  var_pois <- c()
  for (r in 1:nTrials) {
    pois_sample <- rpois(length(LOHcounts_in$n_LOH[idx2]), 
                         mean(LOHcounts_in$n_LOH[idx2]))
    var_pois[r] <- var(pois_sample)
  }
  var_list[[t]] <- var_pois
  p_vals <- rbind(p_vals, c(Tx_name = tx, 
                            p = signif(sum(var_list[[t]] > 
                                             var_sample[[t]])/nTrials, 3)))
  # p_vals <- c(p_vals, sum(var_list[[t]] > var_sample[[t]])/nTrials)
  # names(p_vals[[t]]) <- tx
}

colnames(p_vals) <- c("Tx_name", "p")
p_vals$Tx_name <- factor(p_vals$Tx_name, levels = c("WT", "Cas9", "Drive"))
p_vals$varLabels <- c(paste0(ifelse(p_vals$p != 0,  "var p = ",  "var p < "), 
                             ifelse(p_vals$p != 0, p_vals$p, 1/nTrials)))

# p_vals$Tx_name <- factor(p_vals$Tx_name, levels = c("WT", "Cas9", "Drive"))
Chi_pois_tx <- list()
for (t in 1:3) {
  # t=3
  tx <- levels(nRuns_freq_Tx$Tx_name)[t]
  idx <- nRuns_freq_Tx$Tx_name == tx
  idx2 <- LOHcounts_in$Tx_name == tx
  
  pois_dist <- dpois(nRuns_freq_Tx$value[idx], mean(LOHcounts_in$n_LOH[idx2]), log = F)
  resid <- 1-sum(pois_dist)
  resid_prop <- (pois_dist+resid/length(pois_dist))*resid
  pois_dist_crct <- pois_dist + resid_prop
  Chi_pois_tx[[t]] <- chisq.test(x=nRuns_freq_Tx$cnt[idx], p=pois_dist_crct)
}
Chi_pois_tx[1:3]$p.value
p_vals$chiLabels <- c(paste0("p = ", 
                             signif(c(Chi_pois_tx[[1]]$p.value, 
                                     Chi_pois_tx[[2]]$p.value, 
                                     Chi_pois_tx[[3]]$p.value), 3)))
p_vals$Tx_ID <- substr(p_vals$Tx_name, 1, 1)
p_vals$Tx_ID <- factor(p_vals$Tx_ID, levels = c("D", "C", "W"))

# p_vals$pLabels <- paste0()
nRuns_freq_Tx$Tx_ID <- substr(nRuns_freq_Tx$Tx_name, 1, 1)
nRuns_freq_Tx$Tx_ID <- factor(nRuns_freq_Tx$Tx_ID, levels = c("D", "C", "W"))

LOHcount_isPois <- nRuns_freq_Tx %>% 
  ggplot() + geom_col(aes(x=value, y=prop)) +
  geom_line(data = pois_df, aes(x=value, y=prop), color="red", inherit.aes = F) +
  geom_label(data = p_vals, x = data_len * 0.85, y = 0.15, label.size = 0,
             size = 4, label.padding = unit(0.5, "lines"),
             aes(label = paste0(chiLabels))) +
  xlab("Number of LOH events") + ylab("Proportion") +
  theme(text = element_text(size = 18),
        # strip.text.y = element_text(size = 12),
        # axis.title.x = element_text(size = 12),
        # axis.title.y = element_text(size = 12)
        ) +
  facet_grid(Tx_ID~.)

LOHcount_isPois

ggsave(file.path(outIntDir, "LOHcount_isPois_EC_2022_03.png"), 
       plot = LOHcount_isPois,
       device = "png",
       width = 11, height = 8.5, 
       units = "in",
       dpi = 600)

# Uses old method of estimating n_LOH events without adjusting for complex events ------

# Permutation test of equal mean num of LOH for WT vs ea Tx --------
N_H_perm_list <- all_nGTruns_cln %>% 
  filter(Rep != "00") %>%
  perm_test(df_in = ., cat_var ="Tx_name", cat_names = c("Cas9", "WT"), response_var = "Tot_hom", 
            n_perms = 10000, alpha=0.05, alt_hyp = "two-tailed")

N_F_perm_list <- all_nGTruns_cln %>% 
  filter(Rep != "00") %>%
  perm_test(df_in = ., cat_var ="Tx_name", cat_names = c("Drive", "WT"), response_var = "Tot_hom", 
            n_perms = 10000, alpha=0.05, alt_hyp = "two-tailed")

H_F_perm_list <- all_nGTruns_cln %>% 
  filter(Rep != "00") %>%
  perm_test(df_in = ., cat_var ="Tx_name", cat_names = c("Cas9", "Drive"), response_var = "Tot_hom", 
            n_perms = 10000, alpha=0.05, alt_hyp = "two-tailed")

N_F_perm_list[1:4]

ggplot() + geom_histogram(aes(x=abs(s$permTable[,2]))) + 
  geom_vline(aes(xintercept=s$critVal), color="red") +
  #geom_vline(aes(xintercept=-s$critVal), color="red") +
  geom_vline(aes(xintercept=abs(s$obsvStat)), color="blue")

ggplot(s[[5]], aes(Value)) + stat_ecdf(geom = "step", pad = FALSE)

N_H_test_label <- paste0("mean p = ", round(N_H_perm_list$pVal, 3), 
                         "\nvariance p = ", round(N_H_permVar_list$pVal, 3))

N_F_test_label <- paste0("mean p = ", round(N_F_perm_list$pVal, 3), 
                         "\nvariance p = ", round(N_F_permVar_list$pVal, 3))

H_F_test_label <- paste0("mean p = ", round(H_F_perm_list$pVal, 3), 
                         "\nvariance p = ", round(H_F_permVar_list$pVal, 3))

Tx_LOH_means <- nGTruns_evo_cln %>% 
  dplyr::group_by(Tx_name) %>% 
  dplyr::summarize(mean_LOH = mean(Tot_hom))

lab_y <- max(nGTruns_evo_cln$Tot_hom)

n_LOH_dotplot <- all_nGTruns_cln %>% 
  #filter(!(ID %in% c("N_E09", "N_G03"))) %>%
  filter(Rep != "00") %>% #arrange(desc(Tot_hom))
  #nGTruns_evo_cln %>%
  ggplot() + 
  geom_dotplot(aes(y=Tot_hom, x=Tx_name), binwidth=1, binaxis = "y", 
               stackdir = "center", binpositions = "all",
               dotsize=0.4, stackratio=1) + 
  geom_bracket(xmin = c("WT", "WT", "Cas9"), xmax = c("Cas9", "Drive", "Drive"),
               y.position = c(lab_y + 3, lab_y + 8, lab_y + 13),
               label = c(N_H_test_label, N_F_test_label, H_F_test_label),
               tip.length = 0.025, label.size = 4) +
  stat_summary(aes(y=Tot_hom, x=Tx_name), fun = "mean", fun.min = "mean", fun.max= "mean", 
               size= 0.4, width=0.5, geom = "crossbar") +
  theme_minimal(base_size=14) + 
  theme(panel.grid.minor.x=element_blank(), axis.text.x = ) +
  ylim(c(0, lab_y + 15)) +
  xlab("Drive Type") + ylab("Number of LOH events") 

n_LOH_dotplot

ggsave(file.path(outIntDir, "all_n_LOH_dotplot.png"), 
       plot = n_LOH_dotplot,
       device = "png",
       width = 10, height = 10, 
       units = "in",
       dpi = 300)


# Permutation test of equal variance across Tx -------------

N_H_permVar_list <- all_nGTruns_cln %>% 
  filter(Rep != "00") %>%
  perm_test(data=., cat_var ="Tx_name", cat_names = c("Cas9", "WT"), response_var = "Tot_hom", 
            n_perms = 10000, alpha=0.05, alt_hyp = "two-tailed", test_stat = var)

N_F_permVar_list <- all_nGTruns_cln %>% 
  filter(Rep != "00") %>%
  perm_test(data=., cat_var ="Tx_name", cat_names = c("Drive", "WT"), response_var = "Tot_hom", 
            n_perms = 10000, alpha=0.05, alt_hyp = "two-tailed", test_stat = var)

H_F_permVar_list <- all_nGTruns_cln %>% 
  filter(Rep != "00") %>%
  perm_test(data=., cat_var ="Tx_name", cat_names = c("Cas9", "Drive"), response_var = "Tot_hom", 
            n_perms = 10000, alpha=0.05, alt_hyp = "two-tailed", test_stat = var)


N_H_permVar_list[1:4]
N_F_permVar_list[1:4]
H_F_permVar_list[1:4]

# Test whether LOH accumulation is from Poisson ----------

nRuns_freq <- nGTruns_evo_cln %>% dplyr::count(Tot_hom)
mean_nRuns <- mean(nGTruns_evo_cln$Tot_hom)
pois_dist <- dpois(nRuns_freq$Tot_hom, mean_nRuns, log = F)
resid <- 1-sum(pois_dist)
resid_prop <- (pois_dist+resid/length(pois_dist))*resid
pois_dist_crct <- pois_dist + resid_prop

chi_pois_test <- chisq.test(nRuns_freq$n, p=pois_dist_crct)

ggplot() + geom_col(aes(x=nRuns_freq$Tot_hom, y=nRuns_freq$n/sum(nRuns_freq$n))) +
  geom_line(aes(x=nRuns_freq$Tot_hom, y=pois_dist_crct), color="red")

# 
nRuns_freq_Tx <- nGTruns_evo_cln %>% dplyr::count(Tx_name, Tot_hom)
colnames(nRuns_freq_Tx) <- c("Tx_name", "value", "cnt")
nRuns_freq_Tx <- nRuns_freq_Tx %>% dplyr::group_by(Tx_name) %>% mutate(prop = cnt/sum(cnt))

Chi_pois_tx <- list()
for (t in 1:3) {
  # t=3
  tx <- levels(nRuns_freq_Tx$Tx_name)[t]
  idx <- nRuns_freq_Tx$Tx_name == tx
  idx2 <- nGTruns_evo_cln$Tx_name == tx
  
  pois_dist <- dpois(nRuns_freq_Tx$value[idx], mean(nGTruns_evo_cln$Tot_hom3[idx2]), log = F)
  resid <- 1-sum(pois_dist)
  resid_prop <- (pois_dist+resid/length(pois_dist))*resid
  pois_dist_crct <- pois_dist + resid_prop
  Chi_pois_tx[[t]] <- chisq.test(x=nRuns_freq_Tx$cnt[idx], p=pois_dist_crct)
}

nTrials <- 100000
var_list <- list()
var_sample <- list()
p_vals <- data.frame(NULL)
for (t in 1:3) {
  tx <- levels(nGTruns_evo_cln$Tx_name)[t]
  idx2 <- nGTruns_evo_cln$Tx_name == tx
  var_sample[t] <- var(nGTruns_evo_cln$Tot_hom[idx2])
  var_pois <- c()
  for (r in 1:nTrials) {
    pois_sample <- rpois(length(nGTruns_evo_cln$Tot_hom[idx2]), 
                         mean(nGTruns_evo_cln$Tot_hom[idx2]))
    var_pois[r] <- var(pois_sample)
  }
  var_list[[t]] <- var_pois
  p_vals <- rbind(p_vals, c(Tx_name = tx, p = sum(var_list[[t]] > var_sample[[t]])/nTrials))
  # p_vals <- c(p_vals, sum(var_list[[t]] > var_sample[[t]])/nTrials)
  # names(p_vals[[t]]) <- tx
}

colnames(p_vals) <- c("Tx_name", "p")
p_vals$pLabels <- c(paste0("p < ", 1/nTrials), paste0("p = ", p_vals$p[2:3]))

# Generate Poisson dist for plotting 
pois_dist_tx_crct <- list()
for (t in 1:3) {
  # t = 1
  tx <- levels(nRuns_freq_Tx$Tx_name)[t]
  idx <- nRuns_freq_Tx$Tx_name == tx
  idx2 <- nGTruns_evo_cln$Tx_name == tx
  pois_dist_tx <- dpois(0:max(nRuns_freq_Tx$value), mean(nGTruns_evo_cln$Tot_hom[idx2]), log = F)
  resid <- 1-sum(pois_dist_tx)
  resid_prop <- (pois_dist_tx+resid/length(pois_dist_tx))*resid
  pois_dist_tx_crct[[t]] <- pois_dist_tx + resid_prop
}

pois_df <- data.frame(Tx_name = c(rep("WT", 30), rep("Cas9", 30), rep("Drive", 30)), 
                      value = rep(0:29, 3), prop = unlist(pois_dist_tx_crct))

pois_df$Tx_name <- factor(pois_df$Tx_name, levels = c("WT", "Cas9", "Drive"))
p_vals$Tx_name <- factor(p_vals$Tx_name, levels = c("WT", "Cas9", "Drive"))

nLOH_isPois <- ggplot(data = nRuns_freq_Tx) + geom_col(aes(x=value, y=prop)) +
  geom_line(data = pois_df, aes(x=value, y=prop), color="red", inherit.aes = F) +
  geom_label(data = p_vals, x = 27.5, y = 0.1125, label.size = 0, size = 4,
             aes(label = pLabels)) +
  theme(strip.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12)) +
  #annotate(geom = "text", x = 25, y = 0.10, label=p_vals$p) +
  xlab("# LOH events") + ylab("Proportion") +
  facet_grid(Tx_name~.)

nLOH_isPois

ggsave(file.path(outIntDir, "nLOH_isPois.png"), 
       plot = nLOH_isPois,
       device = "png",
       width = 10, height = 10, 
       units = "in",
       dpi = 300)

