# LOH allele bias deviation from 0.5
#
#
#

# Starts with the dataframes "nGTruns_evo_cln" and "GTruns_evo_cln" generated in <filename>
# which are products of the GTnRunsFunct and GTrunsFunct functions, respectively, 
# followed by filtering to exclude error-prone clones and founders. 

#fract_BY_in <-nGTruns_evo_cln
# fract_BY_in <- all_LOHcounts_GQ70 
# fract_BY_in <- all_longLOHcounts_GQ70 
# all_LOHcounts_err <- all_LOHbounds_err %>% CountLOHevents()



all_LOHbounds_merge_NS %>%
  ggplot() + geom_histogram(aes(x = log10(length)), binwidth = 0.3) + 
  facet_wrap(~GT, ncol = 1, scales = "free_y") +
  xlab("Log(# of markers supporting tract)")

all_LOHbounds_merge %>% 
  filter(length == 1, GT != "0/1", !is_error) %>% 
  count(ID, GT) %>% 
  # arrange(desc(n))
  ggplot() + geom_histogram(aes(x = n, fill = GT), binwidth = 1)
  
fract_BY_in <- all_LOHcounts_merge_NS

fract_BY_in <- all_LOHbounds_merge %>% 
  group_by(Tx_name, ID) %>%
  filter(length == 1) %>% summarize(n_BYtot = sum(GT == "0/0"),
                                    n_RMtot = sum(GT == "1/1"),
                                    n_LOHtot = n())

fract_BY_in$Tx_ID <- Recode_Tx_ID(fract_BY_in$Tx)

# Test whether RM:BY conversion ratio is binomial with Pr(BY)=0.5 ##############
fract_BY_in$n_BYtot <- (fract_BY_in$n_BYsmpl + fract_BY_in$n_BYcmplx)
fract_BY_in$n_RMtot <- (fract_BY_in$n_RMsmpl + fract_BY_in$n_RMcmplx)
fract_BY_in$n_LOHtot <- fract_BY_in$n_BYtot + fract_BY_in$n_RMtot 
# fract_BY_in$n_LOHtot <- (fract_BY_in$n_LOH + fract_BY_in$n_cmplx - fract_BY_in$n_Mix)
# fract_BY_in$n_LOHtot <- (fract_BY_in$n_LOH - fract_BY_in$n_Mix)

# proportion of LOH that are BY allele as an estimate of the Pr = X in Binomial Dist
fract_BY_in$fractLOH_BY <- fract_BY_in$n_BYtot/fract_BY_in$n_LOHtot

mean_fract_BY <- fract_BY_in %>% 
  # group_by(Tx_ID) %>% 
  summarise(grand_mean_f_BY = sum(n_BYtot)/sum(n_LOHtot),
            mean_f_BY = mean(n_BYtot/n_LOHtot),
            f_BY_CIlo = mean(n_BYtot/n_LOHtot) - se(n_BYtot/n_LOHtot) * 1.96,
            f_BY_CIup = mean(n_BYtot/n_LOHtot) + se(n_BYtot/n_LOHtot) * 1.96)

fract_BY_in %>% ggplot() + 
  geom_histogram(aes(x = fractLOH_BY), binwidth = 0.05) + 
  facet_grid(Tx_name~.)

overall_mean <- sum(fract_BY_in$n_BYtot)/sum(fract_BY_in$n_LOHtot)
overall_sd <- sd(fract_BY_in$n_BYtot/fract_BY_in$n_LOHtot, na.rm = T)

LOH_bias_binom <- binom.test(sum(fract_BY_in$n_BYtot), sum(fract_BY_in$n_LOHtot), p = 0.5)

fract_conv_BY <- fract_BY_in %>% 
  dplyr::group_by(Line) %>% 
  dplyr::summarise(n_BY = sum(n_BYtot),
                   n_RM = sum(n_RMtot),
                   n_Tot = sum(n_BYtot) + sum(n_RMtot))

fract_conv_BY %>% ggplot() + geom_point(aes(x = n_BY, y = n_RM)) + geom_abline() + ylim(0, NA) + xlim(0, NA)

binom_df <- data.frame(NULL)
for(l in levels(fract_conv_BY$Line)) {
  # l = "N_A"
  test_row <- fract_conv_BY %>% filter(Line == l)
  binom_result <- binom.test(test_row$n_BY, test_row$n_Tot, alternative = "two.sided")
  binom_p <- binom_result$p.value
  binom_interval <- binom_result$conf.int
  binom_result <- data.frame(Line = l, CI_low = binom_result$conf.int[1],
                             CI_hi = binom_result$conf.int[2],
                             p_val = binom_result$p.value)
  binom_df <- rbind(binom_df, binom_result)
}



fract_conv_BY$fract_BY <- fract_conv_BY$n_BY/fract_conv_BY$n_Tot
mean(fract_conv_BY$fract_BY)
sd(fract_conv_BY$fract_BY)

error_ID_POSi <- all_LOHbounds_merge_EC_50 %>% filter(is_error) %>% 
  paste0(.$ID, "_", .$start_POSi)
all_LOHbounds_50$ID_POSi <- paste0(all_LOHbounds_50$ID, "_", all_LOHbounds_50$POSi)
SNP_bias <- all_LOHbounds_merge_EC_50 %>% filter(GT != "0/1") %>% 
  group_by(GT) %>% summarise(n_SNP = sum(length))
SNP_region_bias <- all_LOHbounds_50 %>% filter(GT != "0/1") %>% 
  group_by(GT) %>% summarize(n_regions = n(), mean_n_SNPs = mean(length))
SNP_merge_bias <- all_LOHbounds_merge_EC_50 %>% filter(GT != "0/1") %>% 
  group_by(GT) %>% summarize(n_regions = n(), mean_n_SNPs = mean(length))

binom.test(x = SNP_bias$n_SNP)

# A simulation to determine whether the data could have come from a binomial distribution
# A number of samples are drawn equal to the number of LOHs in a clone with a Pr(BY) = 0.5
# For each clone, this is repeated nTrials times, and this whole operation is repeated for 
# each clone, which are grouped by Tx. "p_val" returns a list of p-values for each Tx entry.
# For each Tx, the number of BY conversions and the total conversions are each summed, and the 
# estimated Pr(BY) is calculated as sum(BY)/sum(Tot). Likewise, the Pr(BY) for each simulated 
# data set is calculated the same way.

nTrials <- 10000
binom_list <- list()
mean_sample <- list()
p_val <- c()
binom_df <- data.frame(NULL)
for (t in 1:3) {
  # t = 1
  tx <- levels(fract_BY_in$Tx_name)[t]
  nRuns_tx <- fract_BY_in %>% filter(Tx_name == tx)
  n_clns <- nrow(nRuns_tx)
  n_BY <- sum(nRuns_tx$n_BYtot)
  n_conv <- sum(nRuns_tx$n_LOHtot)
  # sum(fract_BY_in$BY[idx2])/sum(fract_BY_in$Tot_hom[idx2])
  mean_sample[[t]] <- mean(nRuns_tx$fractLOH_BY)
  binom_matrix <- matrix(0, nrow = n_clns, ncol = nTrials)
  for (s in 1:n_clns) {
    # s = 1
    binom_matrix[s,] <- rbinom(nTrials, nRuns_tx$n_LOHtot[s], prob = 0.5)
  }
  binom_vals <- apply(binom_matrix, sum, MARGIN = 2)/n_conv
  tx_binom_df <- data.frame(Tx_name = tx, prop_BY = binom_vals)
  binom_df <- rbind(binom_df, tx_binom_df)
  p_val[t] <- sum(binom_vals > n_BY/n_conv)/nTrials
}
binom_df$Tx_name <- factor(binom_df$Tx_name, levels = c("WT", "Cas9", "Drive"))


ggplot() + geom_density(aes(x=binom_vals)) + 
  geom_vline(aes(xintercept = sum(fract_BY_in$n_BYtot)/sum(fract_BY_in$n_LOHtot)), color = "red2") + 
  xlim(c(0, 1))

fract_BY_in %>% ungroup %>% summarise(prop_BY = sum(n_BYtot)/sum(n_LOHtot))

mean_prop_BY <- fract_BY_in %>% 
  group_by(Tx_ID) %>% 
  summarise(grand_mean_f_BY = sum(n_BYtot)/sum(n_LOHtot),
    mean_f_BY = mean(n_BYtot/n_LOHtot),
            f_BY_CIlo = mean(n_BYtot/n_LOHtot) - se(n_BYtot/n_LOHtot) * 1.96,
            f_BY_CIup = mean(n_BYtot/n_LOHtot) + se(n_BYtot/n_LOHtot) * 1.96)

n_c <- max(n_clones_xTx$n)
bi_var <- n_c * 0.25
c_bi <- rbinom(100, n_c, 0.5)
c_bi_var <- var(c_bi)

f_bi_var <- bi_var/n_c^2
f_c_bi <- rbinom(100, n_c, 0.5)/100
f_c_bi_var <- var(f_c_bi)
f_c_bi_var

p_label <- paste0("Binomial test for p = 0.5", paste0("\np = ", p_val))

hist_convBias <- fract_BY_in %>% ggplot() + 
  geom_vline(aes(xintercept = 0.5), size = 0.25) +
  geom_rect(data = mean_prop_BY, aes(xmin = f_BY_CIlo, xmax = f_BY_CIup, 
                                     ymin = -Inf, ymax = Inf), fill = "grey70", alpha = 0.3) +
  geom_histogram(aes(x = fractLOH_BY), binwidth = 0.05) +
  geom_vline(data = mean_prop_BY, aes(xintercept = mean_f_BY), color = "grey10") +
  # geom_label(data = mean_prop_BY, aes(x = mean_f_BY, y = 18, label = signif(mean_f_BY, 3)),
  # label.size = 0, hjust = -0.1, size = 8) +
  # geom_density(data = binom_df, aes(x=prop_BY), color = "red") +
  xlab("Fraction of LOH Conversions to BY Allele") + 
  ylab("Number of Clones") +
  # ylim(0, 20) +
  # xlim(c(0,1)) +
  facet_grid(Tx_ID~.) +
  theme(text = element_text(size = 28),
        panel.grid.minor = element_blank())

hist_convBias

ggsave(file.path(outIntDir, "hist_convBias_2022_03.png"), 
       plot = hist_convBias,
       device = "png",
       width = 10, height = 10, 
       units = "in",
       dpi = 600)

point_convBias <- fract_BY_in %>% 
  ggplot() + 
  geom_jitter(aes(x=n_BYtot, y=n_RMtot), height=0.25, width=0.25, alpha=0.7, size = 0.75) + 
  geom_abline(size = 0.25) +
  xlab("# of Conversions to BY Allele") + ylab("# of Conversions to RM Allele") +
  ylim(c(-1,12)) + xlim(c(-1,12)) +
  theme(panel.background = element_rect(color="grey80"),
        strip.text.x = element_text(size = 12)) +
  facet_grid(.~Tx_name) #+ ggtitle("Distribution of LOH Conversion Bias")

point_convBias

density_convBias <- fract_BY_in %>% 
  ggplot(aes(x=n_BYtot, y=n_RMtot)) + 
  geom_density2d_filled(contour_var = "ndensity", alpha = 0.7, show.legend = F) + 
  scale_fill_manual(values = c("white", "grey90", "grey80", "grey70", 
                               "grey60", "grey50", "grey40", "grey30", 
                               "grey20", "grey10")) +
  geom_jitter(height=0.2, width=0.2, alpha=0.7, size = 0.5, color = "grey10", show.legend = F) + 
  geom_abline(size = 0.25, color = "black") +
  # geom_abline(slope = 1/(1-), size = 0.25, color = "black") +
  # 
  # xlab("# of Conversions to BY Allele") + ylab("# of Conversions to RM Allele") +
  #ylim(c(0,20)) + xlim(c(0,20)) +
  theme(panel.background = element_rect(color="grey80"),
        strip.text.x = element_text(size = 12)) +
  facet_grid(.~Tx_name) #+ ggtitle("Distribution of LOH Conversion Bias")

density_convBias

# Calculate the proportion of conversions to BY for each chromosome and each Tx

chrom_nRuns <- all_LOHbounds_merge %>% 
  ungroup() %>% 
  filter(GT != "0/1", !is_error) %>%
  group_by(Tx_name, CHROM) %>%
  count(GT)

chrom_nRuns$GT <- factor(chrom_nRuns$GT, labels = c("Ref_hom", "Alt_hom"))

chrom_nRuns_wide <- chrom_nRuns %>% 
  pivot_wider(names_from = GT, values_from = n) %>% 
  mutate(Tot = (Ref_hom+Alt_hom), fract_BY = Ref_hom/(Ref_hom+Alt_hom)) %>% ungroup

# Calculate error by assuming data is drawn from a binomial with expected Pr = BY/Tot
bn_t <- function(a, b) {
  binom.test(a, b, p = 0.5, alternative = "two.sided", conf.level = 0.95)$conf_int
}

binom.test(chrom_nRuns_wide$Ref_hom[1], chrom_nRuns_wide$Tot[1])$conf_int

chrom_nRuns_wide <- chrom_nRuns_wide %>% rowwise(Ref_hom, Tot) %>% mutate(bn_lo = binom.test(Ref_hom, Tot)$conf.int[1],
                            bn_hi = binom.test(Ref_hom, Tot)$conf.int[2])

chrom_nRuns_wide$sd <- sqrt(chrom_nRuns_wide$fract_BY*(1-chrom_nRuns_wide$fract_BY)/
                              chrom_nRuns_wide$Tot)
chrom_nRuns_wide$sd_lw <- chrom_nRuns_wide$fract_BY - chrom_nRuns_wide$sd
chrom_nRuns_wide$sd_hi <- chrom_nRuns_wide$fract_BY + chrom_nRuns_wide$sd

plot_bars <- data.frame(from = seq(1.5, 15.5, 2), to = seq(2.5, 16.5, 2))

alleleBias_chrom_point <- ggplot(data = chrom_nRuns_wide, 
                                 aes(x = as.numeric(CHROM) + (as.numeric(Tx_name) - 2)/4, 
                                     y = fract_BY - 0.5)) + 
  annotate(geom = "rect", xmin = plot_bars$from,
           xmax = plot_bars$to, ymin=-Inf, ymax=Inf, alpha=0.2, fill = "grey70") +
  geom_hline(aes(yintercept = 0), size = 0.15, color = "grey50") +
  geom_errorbar(aes(ymin = bn_lo - 0.5, ymax = bn_hi - 0.5), size = 0.3, width = 0.15) +
  geom_point(aes(color = Tx_name), size = 2.5) + 
  geom_text(aes(y = -0.47, label = Tot, color = Tx_name), size = 4) +
  #geom_line(aes(color = Tx_name, group = Tx_name), size = 0.25) + 
  scale_color_manual(values = c("grey20", "gold3", "purple3")) +
  scale_x_continuous(breaks = 1:16) +
  scale_y_continuous(breaks = c(-0.5, 0, 0.5), labels = c("RM", "0", "BY"), limits = c(-0.5, 0.5)) +
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.text.x = element_text(vjust = 10, size = 10),
        axis.ticks.x = element_blank(),
        legend.position = "top") +
  xlab("Chromosome") + ylab("LOH Allele Bias") #+
# facet_grid(Tx_name~.)

alleleBias_chrom_point

ggsave(file.path(outIntDir, "alleleBias_chrom.png"), 
       plot = alleleBias_chrom_point,
       device = "png",
       width = 11, height = 8.5, 
       units = "in",
       dpi = 600)


alleleBias_chrom_point <- ggplot(data = chrom_nRuns_wide, aes(x = chrom_n, y = fract_BY)) + 
  geom_point(aes(color = Tx_name)) + 
  #geom_line(aes(color = Tx_name, group = Tx_name)) + 
  scale_color_manual(values = c("purple3", "gold3",  "grey20")) +
  ylim(c(0,1)) +
  xlab("Chromosome") + ylab("Fraction of LOH with BY Allele") + facet_grid(Tx_name~.)

alleleBias_chrom_point

ggsave(file.path(outIntDir, "alleleBias_chrom.png"), 
       plot = alleleBias_chrom,
       device = "png",
       width = 11, height = 8.5, 
       units = "in",
       dpi = 300)

ggplot(data = chrom_nRuns_wide, aes(x = chrom_n, y = Tx_name)) + 
  geom_tile(aes(height = 1, width = 1, fill = Tot/chrom_lengths)) + 
  geom_label(aes(label = signif(Tot/chrom_lengths, 3)), label.size = 0, alpha = 0.5) #+
# scale_fill_gradient2(high = "orange2", mid = "white", low = "blue3", 
#                      midpoint = 0.5, limits = c(0, 1))
mean(chrom_nRuns_wide$Tot)

# point_convBias <- fract_BY_in %>% 
#   filter(Tot_hom > 0) %>%
#   ggplot() + 
#   geom_jitter(aes(y=BY, x=Tot_hom), height=0.15, width=0.15) + 
#   geom_abline(aes(slope=0.5, intercept=0)) +
#   xlab("Total # of Conversions") + ylab("# of Conversions to BY") +
#   theme(panel.background = element_rect(color="grey80")) +
#   #ylim(c(-1,8)) + xlim(c(-1,16)) +
#   facet_grid(.~Tx_name) + ggtitle("Distribution of LOH Conversions to BY")

ggsave(file.path(outIntDir, "point_convBias.png"), 
       plot = point_convBias,
       device = "png",
       width = 11, height = 8.5, 
       units = "in",
       dpi = 300)

# Permutation test for difference in allele bias among Tx ---------

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
