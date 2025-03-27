

SNPs_merge_finalGT %>%
  filter(ID == "L000") %>%
  ggplot(aes(Ref_DP_final, Alt_DP_final, color = GT)) +
  geom_point(alpha = 0.6) +
  geom_abline() +
  scale_color_manual(values = allelePal, 
                     name = "Genotype", 
                     labels = c("P1 Hom", "Het", "P2 Hom")) +
  labs(x = "# of reads supporting Parent 1 alllele",
       y = "# of reads supporting Parent 2 alllele",
       title = "Callset from genotype reconciliation")

.alpha_H <- 0.05

pars <- data.frame(pi = c(0.03, 0.96, 0.01), mu = c(0.33, 0.5, 0.66), sigma = c(0.1, 0.07, 0.1))

Hets_f_Alt <- SNPs_merge_finalGT %>% 
  filter(QUAL_P1call >= 1000, QUAL_P2call >= 1000, GT == "0/1") %>% 
  pull(f_Alt) %>% 
  DataHist(., seq(0, 1, 0.01))

norm_mix_Hets <- mix(Hets_f_Alt, pars, dist = "norm", 
                     constr = mixconstr(conmu = "MFX", fixmu = c(F, T, F)))

norm_mix_Hets

n_sites_Hets <- SNPs_merge_finalGT %>% 
  filter(QUAL_P1call >= 1000, QUAL_P2call >= 1000, GT == "0/1") %>% 
  nrow()

d1 <-  NormHist(mu = norm_mix_Hets$parameters$mu[1], 
                v = norm_mix_Hets$parameters$sigma[1]^2) *
  n_sites_Hets * norm_mix_Hets$parameters$pi[1]
d2 <-  NormHist(mu = norm_mix_Hets$parameters$mu[2], 
                v = norm_mix_Hets$parameters$sigma[2]^2) * 
  n_sites_Hets * norm_mix_Hets$parameters$pi[2]
d3 <-  NormHist(mu = norm_mix_Hets$parameters$mu[3], 
                v = norm_mix_Hets$parameters$sigma[3]^2) * 
  n_sites_Hets * norm_mix_Hets$parameters$pi[3]

norm_model_Hets <- data.frame(f_Alt = seq(0, 1, 0.01), dist1 = d1, dist2 = d2, dist3 = d3)

high_cut_Hets <- qnorm(0.95, 0.5, norm_mix_Hets$parameters$sigma[2])
low_cut_Hets <- qnorm(0.05, 0.5, norm_mix_Hets$parameters$sigma[2])

SNPs_merge_finalGT <- SNPs_merge_finalGT %>%
  mutate(Cut_hmn = (GT == "0/1" & (f_Alt < low_cut_Hets | f_Alt > high_cut_Hets)))


SNPs_merge_finalGT %>%
  filter(QUAL_P1call >= 1000, QUAL_P2call >= 1000, GT == "0/1") %>%
  ggplot() +
  geom_histogram(aes(f_Alt, fill = Cut_H), binwidth = 0.01) +
  geom_line(data = norm_model_Hets, aes(x = f_Alt, y = dist1 + 1), color = "red3") +
  geom_line(data = norm_model_Hets, aes(x = f_Alt, y = dist2 + 1), color = "blue2") +
  geom_line(data = norm_model_Hets, aes(x = f_Alt, y = dist3 + 1), color = "orange2") +
  scale_fill_manual(values = c("grey40", "black")) +
  xlim(-0.01, 1.01)

SNPs_merge_finalGT %>%
  filter(QUAL_P1call >= 1000, QUAL_P2call >= 1000, GT == "0/1") %>%
  ggplot() +
  geom_histogram(aes(f_Alt, fill = Cut_hmn), binwidth = 0.01) +
  scale_fill_manual(values = c("grey40", "black")) +
  xlim(-0.01, 1.01)


SNPs_merge_finalGT %>%
  filter(ID == "L000") %>%
  filter(GT == "0/1") %>%
  ggplot(aes(Ref_DP_final, Alt_DP_final, color = Cut_hmn)) +
  geom_point(alpha = 0.6) +
  geom_abline() +
  scale_color_manual(values = allelePal, 
                     name = "Genotype", 
                     labels = c("P1 Hom", "Het", "P2 Hom")) +
  labs(x = "# of reads supporting Parent 1 alllele",
       y = "# of reads supporting Parent 2 alllele",
       title = "Callset from genotype reconciliation")

SNPs_merge_finalGT %>%
  filter(ID == "L000") %>%
  filter(GT == "0/1") %>%
  ggplot(aes(Ref_DP_final, Alt_DP_final, color = Cut_H)) +
  geom_point(alpha = 0.6) +
  geom_abline() +
  scale_color_manual(values = allelePal, 
                     name = "Genotype", 
                     labels = c("P1 Hom", "Het", "P2 Hom")) +
  labs(x = "# of reads supporting Parent 1 alllele",
       y = "# of reads supporting Parent 2 alllele",
       title = "Callset from genotype reconciliation")

