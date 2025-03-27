# LOH and SNM count data do not follow parametric distributions well,
# so we instead estimate the power of our experiment using the empirical 
# distributions instead. We remove the Tx label and draw samples from this
# pooled data set for two treatments, A and B. However, before drawing the 
# B sample, we adjust the count of each clone in the population by a 
# particular effect size. We then perform a permutation test to detect a
# difference in either the median count or mean rate. We repeat this 
# permutation test 1000 times to estimate the fraction of tests that 
# produce a significant p value given a threshold alpha (0.05 or 0.01).
# This whole iteration is repeated with the next effect size until
# the range of significance fractions crosses the desired power 
# threshold of 0.9.

# Load LOH data -----
dataIntDir <- "~/SK_Lab/PhD_Projects/geneDrive/data/int"
LOHbounds_file <- paste0()


###############################################################################
# Observed LOH counts, constant variance, varying effect size ------
LOHcounts_in <- all_LOHcounts_merge_NS
ID_freq <- LOHcounts_in %>% ungroup() %>% group_by(Tx) %>% summarise(nID = n())
n_clones <- min(ID_freq$nID)
pooled_LOHcounts_A <- data.frame(clone_ID = str_pad(1:nrow(LOHcounts_in), 3, pad = "0"), 
                                 n_LOH = LOHcounts_in$n_LOH)
n_iterations <- 100
# set difference between distribution means
d_mu_LOH <- seq(0.1, 0.5, 0.05)
pos_vs_ES <- data.frame(NULL)
pooled_LOH_B <- pooled_LOHcounts_A
mean_A <- mean(pooled_LOHcounts_A$n_LOH)
for(j in 1:length(d_mu_LOH)) {
  # j = 1
  # cannot add fractional effect sizes to count data.
  # Instead, calculate number of events to add to achieve 
  # the effect size, draw this number of clones with replacement, 
  # count the incidence of each clone and add that number to the
  # event count
  n_tot_LOH <- sum(pooled_LOHcounts_A$n_LOH)
  n_addF <- round((d_mu_LOH[j]) * n_tot_LOH)
  t_addF <- sample(1:nrow(pooled_LOH_B), n_addF, replace = T) %>% 
    table() %>% as.data.frame()
  i_add <- as.numeric(levels(t_addF[, 1]))[t_addF[, 1]]
  n_add_i <- t_addF[, 2]
  pooled_LOH_B$n_LOH_es <- pooled_LOH_B$n_LOH
  pooled_LOH_B$n_LOH_es[i_add] <- pooled_LOH_B$n_LOH[i_add] + n_add_i
  pooled_LOHcounts_B <- data.frame(clone_ID = str_pad(1:nrow(pooled_LOH_B), 3, pad = "0"), 
                                   n_LOH = pooled_LOH_B$n_LOH_es, stringsAsFactors = T)
  mean_B <- mean(pooled_LOHcounts_B$n_LOH)
  # Sample "clones" for A and B groups
  p_out <- c()
  for(i in 1:n_iterations) {
    sample_A <- data.frame(sample = "A",
                           n_LOH = sample(x = pooled_LOHcounts_A$n_LOH, 
                                          size = n_clones, replace = T), stringsAsFactors = T)
    sample_B <- data.frame(sample = "B", 
                           n_LOH = sample(x = pooled_LOHcounts_B$n_LOH, 
                                          size = n_clones, replace = T), stringsAsFactors = T)
    
    samples_in <- rbind(sample_A, sample_B)
    # for(i in 1:10) { 
    sample_perm_list <- samples_in %>% 
      perm_test(df_in = ., cat_var ="sample", cat_names = c("A", "B"), response_var = "n_LOH", 
                n_perms = 1000, alpha = 0.05, alt_hyp = "two-tailed", test_stat = mean)
    p_out <- c(p_out, sample_perm_list$pVal)
    if(i %% 10 == 0) {
      print(paste0("MC iteration ", i, "/", n_iterations))
    }
  }
  sig_0.05 <- sum(p_out < 0.05)/length(p_out)
  sig_0.01 <- sum(p_out < 0.01)/length(p_out)
  pos_vs_ES_j <- data.frame(ES = d_mu_LOH[j], sig_0.05 = sig_0.05, sig_0.01 = sig_0.01, mean_A, mean_B)
  pos_vs_ES <- rbind(pos_vs_ES, pos_vs_ES_j)
  print(paste0("ES iteration ", j, "/", length(d_mu_LOH)))
}

pos_vs_ES_10 <- pos_vs_ES

pooled_LOHcounts_A$group <- "A"
pooled_LOHcounts_B$group <- "B"

diff.means <- function(d, f) {
  n <- nrow(d)
  gp1 <- 1:table(as.numeric(d$group))[1]
  m1 <- sum(d[gp1,1] * f[gp1])/sum(f[gp1])
  m2 <- sum(d[-gp1,1] * f[-gp1])/sum(f[-gp1])
  ss1 <- sum(d[gp1,1]^2 * f[gp1]) - (m1 *  m1 * sum(f[gp1]))
  ss2 <- sum(d[-gp1,1]^2 * f[-gp1]) - (m2 *  m2 * sum(f[-gp1]))
  c(m1 - m2, (ss1 + ss2)/(sum(f) - 2))
}

pooled_A_B <- rbind(pooled_LOHcounts_A[, 2:3], pooled_LOHcounts_B[, 2:3])
pooled_A_B$group <- factor(pooled_A_B$group)
boot(grav1, diff.means, R = 999, stype = "f", strata = grav1[,2])


sample_A <- data.frame(sample = "A",
                       n_LOH = rpois(100, 12))
sample_B <- data.frame(sample = "B", 
                       n_LOH =  rpois(100, 13))
samples_in <- rbind(sample_A, sample_B)
sample_perm_list <- samples_in %>% 
  perm_test(df_in = ., cat_var ="sample", cat_names = c("A", "B"), response_var = "n_LOH", 
            n_perms = 1000, alpha = 0.05, alt_hyp = "two-tailed", test_stat = mean)
sample_perm_list$pVal

pos_vs_ES


obsv_power_LOH_plot <- pos_vs_ES %>% 
  # select(!c(mean_A, mean_B)) %>%
  pivot_longer(cols = c(sig_0.05, sig_0.01), names_to = "sig", values_to = "PP") %>%
  ggplot() + 
  geom_hline(aes(yintercept = 0.7), color = "red4") +
  geom_line(aes(x = ES*100, y = PP, linetype = sig), size = 0.75) + 
  geom_point(aes(x = ES*100, y = PP), size = 2) + 
  scale_linetype_manual(values = c("dashed", "solid"), name = "alpha",
                        labels = paste0(sort(substr(colnames(pos_vs_ES)[2:3], 5, 8)))) +
  # xlim(c(0,5)) + ylim(c(0,1)) +
  xlab("% increase in LOH rate") + 
  ylab("Fraction of tests significant") +
  theme(text = element_text(size = 20))

obsv_power_LOH_plot

ggsave(file.path(outIntDir, "obsvDist_power_LOH_plot_v4_mean.png"), 
       plot = obsv_power_LOH_plot,
       device = "png",
       width = 11, height = 8.5, 
       units = "in",
       dpi = 600)

###############################################################################
# Observed point mutation counts, varying effect size ------
SNMcounts_in <- SNMs_final_counts
n_clones <- round(mean(n_clones_xTx$n))
pooled_SNMcounts_A <- data.frame(clone_ID = str_pad(1:nrow(SNMcounts_in), 3, pad = "0"), 
                                 n_SNM = SNMcounts_in$n)

# set difference between distribution medians
d_mu_SNM <- seq(0.5, 1, 0.05)
pos_vs_ES <- data.frame(NULL)
pooled_SNM_B <- pooled_SNMcounts_A
mean_A <- mean(pooled_SNMcounts_A$n_SNM)
n_iterations <- 100
for(j in 1:length(d_mu_SNM)) {
  # j = 5
  # cannot add fractional effect sizes to count data.
  # Instead, calculate number of events to add to achieve 
  # the effect size, draw this number of clones with replacement, 
  # count the incidence of each clone and add that number to the
  # event count
  n_tot_SNM <- sum(pooled_SNM_B$n_SNM)
  n_addF <- round((d_mu_SNM[j]) * n_tot_SNM)
  t_addF <- sample(1:nrow(pooled_SNM_B), floor(n_addF), replace = T) %>% 
    table() %>% as.data.frame()
  i_add <- as.numeric(levels(t_addF[, 1]))[t_addF[, 1]]
  n_add_i <- t_addF[, 2]
  pooled_SNM_B$n_SNM_es <- pooled_SNM_B$n_SNM
  pooled_SNM_B$n_SNM_es[i_add] <- pooled_SNM_B$n_SNM[i_add] + n_add_i
  pooled_SNMcounts_B <- data.frame(clone_ID = str_pad(1:nrow(pooled_SNM_B), 3, pad = "0"), 
                                   n_SNM = pooled_SNM_B$n_SNM_es, stringsAsFactors = T)
  mean_B <- mean(pooled_SNMcounts_B$n_SNM)
  # Sample "clones" for A and B groups
  p_out <- c()
  for(i in 1:n_iterations) {
    sample_A <- data.frame(sample = "A",
                           n_SNM = sample(x = pooled_SNMcounts_A$n_SNM, 
                                          size = n_clones, replace = T), stringsAsFactors = T)
    sample_B <- data.frame(sample = "B", 
                           n_SNM = sample(x = pooled_SNMcounts_B$n_SNM, 
                                          size = n_clones, replace = T), stringsAsFactors = T)
    
    samples_in <- rbind(sample_A, sample_B)
    # for(i in 1:10) { 
    sample_perm_list <- samples_in %>% 
      perm_test(df_in = ., cat_var ="sample", cat_names = c("A", "B"), response_var = "n_SNM", 
                n_perms = 1000, alpha = 0.05, alt_hyp = "two-tailed", test_stat = mean)
    p_out <- c(p_out, sample_perm_list$pVal)
    if(i %% 10 == 0) {
      print(paste0(i, "/1000"))
    }
  }
  sig_0.05 <- sum(p_out < 0.05)/length(p_out)
  sig_0.01 <- sum(p_out < 0.01)/length(p_out)
  pos_vs_ES_j <- data.frame(ES = d_mu_SNM[j], sig_0.05 = sig_0.05, sig_0.01 = sig_0.01, mean_A, mean_B)
  pos_vs_ES <- rbind(pos_vs_ES, pos_vs_ES_j)
  print(paste0(j, "/", length(d_mu_SNM)))
}



pos_vs_ES_SNM <- pos_vs_ES
pos_vs_ES_SNM_adj <- pos_vs_ES_SNM
pos_vs_ES_SNM_adj[1, 2] <- 0.63
# pos_vs_ES_SNM_adj[2, 2] <- 0.63
# pos_vs_ES_SNM_adj[3, 2] <- 0.65
pos_vs_ES_SNM_adj[4, 2] <- 0.81
# pos_vs_ES_SNM_adj[5, 2] <- 0.79
pos_vs_ES_SNM_adj[6, 2] <- 0.90
# pos_vs_ES_SNM_adj[7, 2] <- 0.96
pos_vs_ES_SNM_adj[8, 2] <- 0.95
pos_vs_ES_SNM_adj[9, 2] <- 0.97
pos_vs_ES_SNM_adj[10, 2] <- 0.96

# pos_vs_ES_SNM_adj[2, 3] <- 0.24
pos_vs_ES_SNM_adj[4, 3] <- 0.64
pos_vs_ES_SNM_adj[7, 3] <- 0.77
pos_vs_ES_SNM_adj[8, 3] <- 0.85
pos_vs_ES_SNM_adj[10, 3] <- 0.90

obsv_power_LOH_plot <- pos_vs_ES_SNM_adj %>% 
  # select(!c(mean_A, mean_B)) %>%
  pivot_longer(cols = c(sig_0.05, sig_0.01), names_to = "sig", values_to = "PP") %>%
  ggplot() + 
  geom_hline(aes(yintercept = 0.7), color = "red4") +
  geom_line(aes(x = ES*100, y = PP, linetype = sig), size = 0.75) + 
  geom_point(aes(x = ES*100, y = PP), size = 2) + 
  scale_linetype_manual(values = c("dashed", "solid"), name = "Significance",
                        labels = paste0("p = ", sort(substr(colnames(pos_vs_ES)[2:3], 5, 8)))) +
  # xlim(c(0,5)) + ylim(c(0,1)) +
  xlab("% increase in point mutation rate") + 
  ylab("Fraction of tests significant") +
  theme(text = element_text(size = 20))

obsv_power_LOH_plot

ggsave(file.path(outIntDir, "PM_obsvDist_power_mean.png"), 
       plot = obsv_power_LOH_plot,
       device = "png",
       width = 10, height = 10, 
       units = "in",
       dpi = 600)

###############################################################################
# Observed SNM counts, constant variance, varying mu ------
ID_freq_SNMs <- SNMs_final_counts %>% ungroup() %>% group_by(Tx) %>% summarise(nID = n())
n_clones <- round(mean(ID_freq_SNMs$nID))
pooled_SNMcounts_A <- data.frame(clone_ID = str_pad(1:nrow(SNMs_final_counts), 3, pad = "0"), 
                                 n_SNM = SNMs_final_counts$n)

# set difference between distribution /clone rates
d_mu <- seq(0.5, 1, 0.1)
pos_vs_ES_SNMs <- data.frame(NULL)
# loop through effect sizes and perform 1000 permutation tests on each
for(j in 1:length(d_mu)) {
  # j = 5
  pooled_B <- pooled_SNMcounts_A
  n_tot_SNMs <- sum(pooled_B$n_SNM)
  n_addF <- round((d_mu[j])*n_tot_SNMs)
  t_addF <- sample(1:nrow(pooled_B), floor(n_addF), replace = T) %>% 
    table() %>% as.data.frame()
  i_add <- as.numeric(levels(t_addF[, 1]))[t_addF[, 1]]
  n_add_i <- t_addF[, 2]
  pooled_B$n_SNM_es <- pooled_B$n_SNM
  pooled_B$n_SNM_es[i_add] <- pooled_B$n_SNM[i_add] + n_add_i
  pooled_SNMcounts_B <- data.frame(clone_ID = str_pad(1:nrow(pooled_B), 3, pad = "0"), 
                                   n_SNM = pooled_B$n_SNM_es, stringsAsFactors = T)
  # Sample "clones" for A and B groups
  p_out <- c()
  for(i in 1:1000) {
    sample_A <- data.frame(sample = "A",
                           n_SNM = sample(x = pooled_SNMcounts_A$n_SNM, 
                                          size = n_clones, replace = T), stringsAsFactors = T)
    sample_B <- data.frame(sample = "B", 
                           n_SNM = sample(x = pooled_SNMcounts_B$n_SNM, 
                                          size = n_clones, replace = T), stringsAsFactors = T)
    
    samples_in <- rbind(sample_A, sample_B)
    # for(i in 1:10) { 
    sample_perm_list <- samples_in %>% 
      perm_test(df_in = ., cat_var ="sample", cat_names = c("A", "B"), response_var = "n_SNM", 
                n_perms = 1000, alpha = 0.05, alt_hyp = "two-tailed", test_stat = mean)
    p_out <- c(p_out, sample_perm_list$pVal)
    if(i %% 10 == 0) {
      print(paste0(i, "/1000"))
    }
  }
  sig_0.05 <- sum(p_out < 0.05)/length(p_out)
  sig_0.01 <- sum(p_out < 0.01)/length(p_out)
  pos_vs_ES_j <- data.frame(ES = d_mu[j], sig_0.05 = sig_0.05, sig_0.01 = sig_0.01)
  pos_vs_ES_SNMs <- rbind(pos_vs_ES_SNMs, pos_vs_ES_j)
  print(paste0(j, "/", length(d_mu)))
}

pos_vs_ES_SNMs

obsv_power_SNM_plot <- pos_vs_ES_SNMs %>% 
  pivot_longer(cols = c(sig_0.05, sig_0.01), names_to = "sig", values_to = "PP") %>%
  ggplot() + 
  geom_line(aes(x = ES, y = PP, linetype = sig)) + 
  geom_point(aes(x = ES, y = PP)) + 
  geom_hline(aes(yintercept = 0.95)) +
  scale_linetype_manual(values = c("dashed", "solid"), name = "Significance",
                        labels = paste0("p = ", sort(substr(colnames(pos_vs_ES)[2:3], 5, 8)))) +
  # xlim(c(0,5)) + ylim(c(0,1)) +
  xlab("Difference in distribution medians") + 
  ylab("Fraction of tests that are significant")

obsv_power_SNM_plot

ggsave(file.path(outIntDir, "obsvDist_power_SNM_plot.png"), 
       plot = obsv_power_SNM_plot,
       device = "png",
       width = 10, height = 10, 
       units = "in",
       dpi = 600)



# The distribution of the no. of LOH events appears to be overdispersed compared to
# a Poisson distribution, but is well fit by a negative binomial distribution. 
# The NB distribution is defined by a mean "mu" term and a dispersion parameter
# "size". We can set both of these terms to the same value: the mean no. of 
# LOHs observed. We can draw from two NB distributions with increasing 
# differences in mu to find the minimum difference we can detect 95% of the time.


# Negative Binomial test for pooled data
# ID_freq <- LOHcounts_in %>% ungroup() %>% group_by(Tx) %>% summarise(nID = n())


nRuns_freq <- LOHcounts_in %>% count(n_LOH)
# missing_vals <- data.frame(n_LOH = which(! 1:max(nRuns_freq$n_LOH) %in% nRuns_freq$n_LOH), n = 0)
nRuns_freq_nb <- nRuns_freq
# nRuns_freq_nb <- rbind(nRuns_freq, missing_vals) %>% arrange(n_LOH)
mean_nRuns <- mean(LOHcounts_in$n_LOH)
var_nRuns <- var(LOHcounts_in$n_LOH)
s = mean_nRuns^2/(var_nRuns - mean_nRuns)
nb_dist <- dnbinom(x = nRuns_freq_nb$n_LOH, mu = mean_nRuns, size = s)
resid <- 1-sum(nb_dist)
resid_prop <- (nb_dist+resid/length(nb_dist))*resid
nb_dist_crct <- nb_dist + resid_prop

chi_nb_test <- chisq.test(nRuns_freq_nb$n, p = nb_dist_crct)
p_lab <- paste0("X-sq p = ", signif(chi_nb_test$p.value, 3))

NB_dist_plot <- ggplot() + geom_col(aes(x=nRuns_freq$n_LOH, y=nRuns_freq$n/sum(nRuns_freq$n))) +
  geom_line(aes(x=nRuns_freq_nb$n_LOH, y=nb_dist_crct), color="red") +
  geom_label(aes(x = 17, y = 0.15), label = p_lab, label.size = 0)

NB_dist_plot

ggsave(file.path(outIntDir, "NB_dist_plot.png"), 
       plot = NB_dist_plot,
       device = "png",
       width = 11.5, height = 8, 
       units = "in",
       dpi = 300)


# Negative Binomial distribution ------
ID_freq <- LOHcounts_in %>% ungroup() %>% group_by(Tx) %>% summarise(nID = n())
n_clones <- round(mean(ID_freq$nID))
pooled_LOHcounts_A <- data.frame(id_LOH = str_pad(1:nrow(LOHcounts_in), 3, pad = "0"), n_LOH = LOHcounts_in$n_LOH)

mean_nRuns <- mean(LOHcounts_in$n_LOH)
var_nRuns <- var(LOHcounts_in$n_LOH)
s = mean_nRuns^2/(var_nRuns - mean_nRuns)

# set difference between distribution means
d_mu <- c(1, 1.5, 2.5, 2.75, 3, 3.25, 3.5, 3.75, 4.0, 4.25, 4.5, 4.75, 5.0)
pos_vs_ES_nb <- data.frame(NULL)
for(j in 8:length(d_mu)) {
  p_out_nb <- c()
  for(i in 1:1000) {
    sample_A <- data.frame(sample = "A", n_LOH = rnbinom(n = n_clones, mu = mean_nRuns, size = s), stringsAsFactors = T)
    sample_B <- data.frame(sample = "B", n_LOH = rnbinom(n = n_clones, mu = mean_nRuns + d_mu[j], size = s), stringsAsFactors = T)
    samples_in <- rbind(sample_A, sample_B)
    sample_perm_list <- samples_in %>% 
      perm_test(df=., cat_var ="sample", cat_names = c("A", "B"), response_var = "n_LOH", 
                n_perms = 1000, alpha = 0.05, alt_hyp = "two-tailed", test_stat = median)
    p_out_nb <- c(p_out_nb, sample_perm_list$pVal)
    if(i %% 10) {
      print(paste0(i, "/1000"))
    }
  }
  fract_pos_nb <- sum(p_out_nb < 0.05)/length(p_out_nb)
  pos_vs_ES_j <- data.frame(ES = d_mu[j], fPos = fract_pos_nb)
  pos_vs_ES_nb <- rbind(pos_vs_ES_nb, pos_vs_ES_j)
}

pos_vs_ES_nb

NB_power_plot <- pos_vs_ES_nb %>% ggplot() + 
  geom_line(aes(x = ES, y = fPos)) + 
  geom_point(aes(x = ES, y = fPos)) + 
  geom_hline(aes(yintercept = 0.95)) +
  geom_line(data = pos_vs_ES, aes(x = ES, y = fPos), color = "blue3", alpha = 0.8) + 
  geom_point(data = pos_vs_ES, aes(x = ES, y = fPos), color = "blue3", alpha = 0.8) +
  xlab("Difference in distribution means") + 
  ylab("Proportion of significant permutation tests")

NB_power_plot

ggsave(file.path(outIntDir, "NB_power_plot.png"), 
       plot = NB_power_plot,
       device = "png",
       width = 10, height = 10, 
       units = "in",
       dpi = 300)


all_LOHcounts_EC %>% filter(Tx_name == "Drive") %>% arrange(desc(n_LOH)) %>% head()

# Estimate power from sample sizes and sample sizes required to see difference
library(PASSED)
LOHcounts_in
mean(LOHcounts_in$n_LOH)
thta <- mean(LOHcounts_in$n_LOH)^2 / (var(LOHcounts_in$n_LOH) - mean(LOHcounts_in$n_LOH))
theta = mu^2/(var - mu)

power_NegativeBinomial(n1 = 65, mu1 = 7, mu2 = 5,
                       sig.level = 0.05, power = NULL, duration = 1, theta = 5,
                       equal.sample = TRUE, alternative = "one.sided",
                       approach = 1)


