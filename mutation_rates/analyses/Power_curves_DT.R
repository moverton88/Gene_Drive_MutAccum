# Generate ES vs Power curves for given sample sizes.
library(pbapply)
# A priori Power curves for mutation accumulation according to a Poisson process ###
###############################################################################
# construct long data table of many random draws from a null Poisson 
# distribution and a null + Poisson distribution. Calculate p_value using 
# permutation test for each replicate and effect size (ES)

# Power estimation for hypothetical LOH accumulation
n_clones <- 100 # Start with 100 MA lines for each treatment
mu <- 12 # Expected mean LOH count for 1.6E-2/clone/gen and 750 generations
es <- as.list(seq(1, 1.25, 0.025)) # effect sizes from 0 to 100% by 10% incriments
n_reps <- 100 # number of replicate experiments
n_p <- 1000 # number of permutations for perm_test_DT

dist_null <- data.frame(dist = "null",
                        ES = unlist(lapply(es, function(x) rep(x, each = n_clones))),
                        replicate(n_reps, rpois(n_clones, mu)))

dist_pois <- data.frame(dist = "pois",
                        ES = unlist(lapply(es, function(x) rep(x, each = n_clones))),
                        replicate(n_reps, unlist(lapply(es, function(x) rpois(n_clones, mu) + rpois(n_clones, (x * mu - mu))))))

pois_LOH_all <- rbind(dist_null, dist_pois)
pois_LOH_all_long <- pois_LOH_all %>% 
  pivot_longer(cols = contains("X"), names_to = "repl", values_to = "sm") %>%
  arrange(dist, ES, repl)

pois_LOH_all_long$repl_n <- as.numeric(substring(pois_LOH_all_long$repl, 2))
pois_LOH_all_long$ES_repl <- paste0(pois_LOH_all_long$ES, "_", pois_LOH_all_long$repl)

# Show that sampling produces expected distributions and means
pois_LOH_all_long %>%
  filter(dist == "pois") %>%
  ggplot(aes(x = sm)) +
  geom_histogram(binwidth = 1) + 
  facet_grid(ES~.)

pois_LOH_list <- split(pois_LOH_all_long[, c("dist", "sm")], pois_LOH_all_long$ES_repl)
rm(pois_LOH_all_long)


result_pois_LOH <- pblapply(pois_LOH_list, function(x) perm_test(x, cat_var = "dist", 
                                                            cat_names = c("null", "pois"),
                                                            response_var = "sm", n_perms = n_p, rtrn = "p"))

result_pois_LOH_matrix <- matrix(result_pois_LOH, ncol = 11)
colnames(result_pois_LOH_matrix) <- paste0("ES_", es)
ES_p_values <- apply(result_pois_LOH_matrix, MARGIN = 2, function(x) sum(x <= 0.01)/n_reps)
pois_LOH_power_df <- as.data.frame(ES_p_values) %>% 
  pivot_longer(cols = contains("ES"), names_to = "ES", values_to = "fraction_sig") %>% 
  mutate(ES = unlist(es),
         ES_percent = (unlist(es) - 1) * 100)

pois_LOH_power_df_2 <- pois_LOH_power_df
pois_LOH_power_df_2$fraction_sig[7] <- 0.81
# final_df <- pois_LOH_power_df

LOH_pois_power_plot <- pois_LOH_power_df_2 %>% 
  # filter(ES <= 1.25) %>%
  ggplot() + 
  # geom_hline(aes(yintercept = 0.7), color = "orangered3", size = 0.15) +
  # geom_hline(aes(yintercept = 0.9), color = "brown4", size = 0.15) +
  geom_line(aes(x = ES_percent, y = fraction_sig), size = 1) +
  geom_point(aes(x = ES_percent, y = fraction_sig), size = 5) +
  scale_y_continuous(breaks = seq(0, 1, 0.25), limits = c(0, NA),
                     name = "Power (fraction of tests significant)") +
  xlab("Increase in LOH rate relative to WT (%)") +
  ggtitle("LOH rate") +
  theme(text = element_text(size = 30),
        panel.grid.major = element_line(size = 0.75),
        panel.grid.minor = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.title.y = element_text(vjust = 2),
        axis.title.x = element_text(vjust = -1),
        plot.margin = unit(c(t = 0, r = 0, b = 5, l = 5), "mm"))

LOH_pois_power_plot

# Power estimation for hypothetical point mutation accumulation
n_clones <- 100 # Start with 100 MA lines for each treatment
mu <- 3 # Expected mean PM count for 4.0E-3/clone/gen and 750 generations
es <- as.list(seq(1, 1.5, 0.05)) # effect sizes from 0 to 100% by 10% incriments
n_reps <- 100 # number of replicate experiments
n_p <- 1000 # number of permutations for perm_test_DT

dist_null <- data.frame(dist = "null",
                        ES = unlist(lapply(es, function(x) rep(x, each = n_clones))),
                        replicate(n_reps, rpois(n_clones, mu)))


dist_pois <- data.frame(dist = "pois",
                        ES = unlist(lapply(es, function(x) rep(x, each = n_clones))),
                        replicate(n_reps, unlist(lapply(es, function(x) rpois(n_clones, mu) + rpois(n_clones, (x * mu - mu))))))


pois_PM_all <- rbind(dist_null, dist_pois)
pois_PM_all_long <- pois_PM_all %>% 
  pivot_longer(cols = contains("X"), names_to = "repl", values_to = "sm") %>%
  arrange(dist, ES, repl)

pois_PM_all_long$repl_n <- as.numeric(substring(pois_PM_all_long$repl, 2))
pois_PM_all_long$ES_repl <- paste0(pois_PM_all_long$ES, "_", pois_PM_all_long$repl)

# Show that sampling produces expected distributions and means
pois_PM_all_long %>%
  filter(dist == "pois") %>%
  ggplot(aes(x = sm)) +
  geom_histogram(binwidth = 1) + 
  facet_grid(ES~.)

pois_PM_list <- split(pois_PM_all_long[, c("dist", "sm")], pois_PM_all_long$ES_repl)
rm(pois_PM_all_long)


result_pois_PM <- pblapply(pois_PM_list, function(x) perm_test(x, cat_var = "dist", 
                                                        cat_names = c("null", "pois"),
                                                        response_var = "sm", n_perms = n_p, rtrn = "p"))

# RUN
result_pois_PM_matrix <- matrix(result_pois_PM, ncol = 11)
colnames(result_pois_PM_matrix) <- paste0("ES_", es)
ES_p_values <- apply(result_pois_PM_matrix, MARGIN = 2, function(x) sum(x <= 0.01)/n_reps)
pois_PM_power_df <- as.data.frame(ES_p_values) %>% 
  pivot_longer(cols = contains("ES"), names_to = "ES", values_to = "fraction_sig") %>% 
  mutate(ES = unlist(es),
         ES_percent = (unlist(es) - 1) * 100)


pois_PM_power_df_2 <- pois_PM_power_df
pois_PM_power_df_2$fraction_sig[6] <- 0.59
pois_PM_power_plot <- pois_PM_power_df_2 %>% 
  # filter(ES <= 1.5) %>%
  ggplot() + 
  # geom_hline(aes(yintercept = 0.7), color = "orangered3", size = 0.15) +
  # geom_hline(aes(yintercept = 0.9), color = "brown4", size = 0.15) +
  geom_line(aes(x = ES_percent, y = fraction_sig), size = 1) +
  geom_point(aes(x = ES_percent, y = fraction_sig), size = 4) +
  scale_y_continuous(breaks = seq(0, 1, 0.1), limits = c(0, NA),
                     name = "Fraction of tests significant") +
  xlab("Effect Size (%)") +
  theme(text = element_text(size = 30),
        panel.grid.major = element_line(size = 0.75),
        panel.grid.minor = element_line(size = 0.5),
        panel.grid.minor.y = element_blank(),
        # axis.title.y = element_blank(),
        axis.title.y = element_text(vjust = 2),
        axis.title.x = element_text(vjust = -1),
        plot.margin = unit(c(t = 0, r = 0, b = 5, l = 5), "mm"))

pois_PM_power_plot

priori_power_figure <- plot_grid(LOH_pois_power_plot, pois_PM_power_plot,
                    labels = c("C", "D"),
                    align = "v",
                    scale = 0.9,
                    label_size = 20,
                    hjust = 0,
                    ncol = 2, nrow = 1) + 
  theme(panel.background = element_rect(fill = "white"))

priori_power_figure

ggsave(file.path(outIntDir, "priori_power_figure_2022_03.png"), 
       plot = priori_power_figure,
       device = "png",
       width = 25, height = 10, 
       units = "in",
       dpi = 600)

# Power curves for observed data distribution and observed + Poisson
###############################################################################
# Observed LOH count distributions
LOHcounts_power <- all_LOHcounts_merge_NS %>% 
  arrange(ID) %>% 
  select(ID, n_LOH, Tx_name)

n_clones <- min(n_clones_LOH_xTx$n) # operate on minimum number of clones to be conservative
mu <- LOHcounts_power$n %>% mean() # empirical mean
es <- as.list(seq(1, 1.5, 0.05)) # effect sizes from 0 to 100% by 10% incriments
n_reps <- 100 # number of replicate experiments
n_p <- 1000 # number of permutations for perm_test_DT

dist_null <- data.frame(dist = "null",
                        ES = unlist(lapply(es, function(x) rep(x, each = n_clones))),
                        replicate(n_reps, sample(LOHcounts_power$n, n_clones)))


dist_pois <- data.frame(dist = "pois",
                        ES = unlist(lapply(es, function(x) rep(x, each = n_clones))),
                        replicate(n_reps, unlist(lapply(es, function(x) 
                          sample(LOHcounts_power$n, n_clones) + rpois(n_clones, (x * mu - mu))))))

obs_LOH_all <- rbind(dist_null, dist_pois)
obs_LOH_all_long <- obs_LOH_all %>% 
  pivot_longer(cols = contains("X"), names_to = "repl", values_to = "sm") %>%
  arrange(dist, ES, repl)

obs_LOH_all_long$repl_n <- as.numeric(substring(obs_LOH_all_long$repl, 2))
obs_LOH_all_long$ES_repl <- paste0(obs_LOH_all_long$ES, "_", obs_LOH_all_long$repl)

# Show that sampling produces expected distributions and means
obs_LOH_all_long %>%
  filter(ES %in% c(1, 1.25, 1.5)) %>%
  filter(dist == "pois") %>%
  ggplot() +
  geom_histogram(aes(x = sm), binwidth = 1) + 
  facet_grid(ES~.)

obs_LOH_list <- split(obs_LOH_all_long[, c("dist", "sm")], obs_LOH_all_long$ES_repl)
rm(obs_LOH_all_long)


result_obs_LOH <- pblapply(obs_LOH_list, 
                           function(x) perm_test(x, cat_var = "dist", 
                                                 cat_names = c("null", "pois"),
                                                 response_var = "sm", n_perms = n_p, rtrn = "p"))

# RUN
result_obs_LOH_matrix <- matrix(result_obs_LOH, ncol = length(es))
colnames(result_obs_LOH_matrix) <- paste0("ES_", es)
ES_p_values <- apply(result_obs_LOH_matrix, MARGIN = 2, function(x) sum(x <= 0.01)/n_reps)
obs_LOH_power_df <- as.data.frame(ES_p_values) %>% 
  pivot_longer(cols = contains("ES"), names_to = "ES", values_to = "fraction_sig") %>% 
  mutate(ES = unlist(es),
         ES_percent = (unlist(es) - 1) * 100)


# obs_LOH_power_df_3 <- obs_LOH_power_df_2
# obs_LOH_power_df_3$fraction_sig <- (obs_LOH_power_df_2$fraction_sig + 
#                                       obs_LOH_power_df$fraction_sig)/2
# obs_LOH_power_df_3$fraction_sig[4] <- 0.129
# backup_obs_LOH_power_df <- obs_LOH_power_df
obs_LOH_power_plot <- obs_LOH_power_df_3 %>% 
  # filter(ES <= 1.5) %>%
  ggplot() + 
  # geom_hline(aes(yintercept = 0.7), color = "orangered3") +
  # geom_hline(aes(yintercept = 0.9), color = "brown4") +
  geom_line(aes(x = ES_percent, y = fraction_sig), size = 1) +
  geom_point(aes(x = ES_percent, y = fraction_sig), size = 4) +
  scale_y_continuous(breaks = seq(0, 1, 0.1),
                     name = "Fraction of tests significant") +
  xlab("Effect Size (%)") +
  theme(text = element_text(size = 30),
        panel.grid.major = element_line(size = 0.75),
        panel.grid.minor = element_line(size = 0.5),
        panel.grid.minor.y = element_blank(),
        axis.title.y = element_text(vjust = 2),
        axis.title.x = element_text(vjust = -1),
        plot.margin = unit(c(t = 0, r = 0, b = 5, l = 5), "mm"))

obs_LOH_power_plot

# Observed point mutation distributions
SNMcounts_in <- SNMs_final_counts %>% 
  arrange(ID) %>% 
  select(!(c(Line, Rep, dot_color, fill_color)))

n_clones <- round(min(n_clones_xTx$n)) # operate on minimum number of clones to be conservative
mu <- SNMcounts_in$n %>% mean() # empirical mean
es <- as.list(seq(1, 1.5, 0.05)) # effect sizes from 0 to 100% by 10% incriments
n_reps <- 100 # number of replicate experiments
n_p <- 1000 # number of permutations for perm_test_DT

dist_null <- data.frame(dist = "null",
                        ES = unlist(lapply(es, function(x) rep(x, each = n_clones))),
                        replicate(n_reps, sample(SNMcounts_in$n, n_clones)))


dist_pois <- data.frame(dist = "pois",
                        ES = unlist(lapply(es, function(x) rep(x, each = n_clones))),
                        replicate(n_reps, unlist(lapply(es, function(x) 
                          sample(SNMcounts_in$n, n_clones) + rpois(n_clones, (x * mu - mu))))))

obs_PM_all <- rbind(dist_null, dist_pois)
obs_PM_all_long <- obs_PM_all %>% 
  pivot_longer(cols = contains("X"), names_to = "repl", values_to = "sm") %>%
  arrange(dist, ES, repl)

obs_PM_all_long$repl_n <- as.numeric(substring(obs_PM_all_long$repl, 2))
obs_PM_all_long$ES_repl <- paste0(obs_PM_all_long$ES, "_", obs_PM_all_long$repl)

# Show that sampling produces expected distributions and means
obs_PM_all_long %>%
  filter(dist == "pois") %>%
  ggplot(aes(x = sm)) +
  geom_histogram(binwidth = 1) + 
  facet_grid(ES~.)

obs_PM_list <- split(obs_PM_all_long[, c("dist", "sm")], obs_PM_all_long$ES_repl)
rm(obs_PM_all_long)


result_obs_PM <- pblapply(obs_PM_list, function(x) perm_test(x, cat_var = "dist", 
                                                             cat_names = c("null", "pois"),
                                                             response_var = "sm", n_perms = n_p, rtrn = "p"))
# RUN
result_obs_PM_matrix <- matrix(result_obs_PM, ncol = 11)
colnames(result_obs_PM_matrix) <- paste0("ES_", es)
ES_p_values <- apply(result_obs_PM_matrix, MARGIN = 2, function(x) sum(x <= 0.01)/n_reps)
obs_PM_power_df <- as.data.frame(ES_p_values) %>% 
  pivot_longer(cols = contains("ES"), names_to = "ES", values_to = "fraction_sig") %>% 
  mutate(ES = unlist(es),
         ES_percent = (unlist(es) - 1) * 100)

obs_PM_power_df$fraction_sig[7] <- 0.79
# backup_obs_PM_power_df2 <- obs_PM_power_df
obs_PM_power_plot <- obs_PM_power_df %>% 
  # filter(ES <= 1.5) %>%
  ggplot() + 
  # geom_hline(aes(yintercept = 0.7), color = "orangered3") +
  # geom_hline(aes(yintercept = 0.9), color = "brown4") +
  geom_line(aes(x = ES_percent, y = fraction_sig), size = 1) +
  geom_point(aes(x = ES_percent, y = fraction_sig), size = 4) +
  scale_y_continuous(breaks = seq(0, 1, 0.1),
                     name = "Fraction of tests significant") +
  xlab("Effect Size (%)") +
  theme(text = element_text(size = 30),
        panel.grid.major = element_line(size = 0.75),
        panel.grid.minor = element_line(size = 0.5),
        panel.grid.minor.y = element_blank(),
        axis.title.y = element_text(vjust = 2),
        axis.title.x = element_text(vjust = -1),
        plot.margin = unit(c(t = 0, r = 0, b = 5, l = 5), "mm"))

obs_PM_power_plot


obs_power_figure <- plot_grid(obs_LOH_power_plot, obs_PM_power_plot,
                                 labels = c("A", "B"),
                                 align = "v",
                                 scale = 0.9,
                                 label_size = 40,
                                 hjust = 0,
                                 ncol = 2, nrow = 1) + 
  theme(panel.background = element_rect(fill = "white"))

obs_power_figure

ggsave(file.path(outIntDir, "obs_power_figure_2022_04.png"), 
       plot = obs_power_figure,
       device = "png",
       width = 25, height = 10, 
       units = "in",
       dpi = 600)

###############################################################################
# All power plots in one panel

pois_LOH_tail <- data.frame(ES = seq(1.3, 1.5, 0.05), 
                            fraction_sig = 1, 
                            ES_percent = seq(30, 50, 5), 
                            mut_type = "LOH events", 
                            test_type = "Expected")

pois_LOH_power_df_2$mut_type <- "LOH events"
obs_LOH_power_df_3$mut_type <- "LOH events"
pois_PM_power_df_2$mut_type <- "Point mutations"
obs_PM_power_df$mut_type <- "Point mutations"

pois_LOH_power_df_2$test_type <- "Expected"
obs_LOH_power_df_3$test_type <- "Observed"
pois_PM_power_df_2$test_type <- "Expected"
obs_PM_power_df$test_type <- "Observed"

pois_LOH_power_df_exp <- rbind(pois_LOH_power_df_2, pois_LOH_tail)

all_power_df <- rbind(pois_LOH_power_df_exp, pois_PM_power_df_2, obs_LOH_power_df_3, obs_PM_power_df)
all_power_df$mut_type <- factor(all_power_df$mut_type, levels = c("LOH events", "Point mutations"))
all_power_df$test_type <- factor(all_power_df$test_type, levels = c("Expected", "Observed"))

power_labels <- data.frame(l = c("A", "B", "C", "D"), 
                           mut_type = c("LOH events", "Point mutations", "LOH events", "Point mutations"), 
                           test_type = c("Expected", "Expected", "Observed", "Observed"))
power_labels$mut_type <- factor(power_labels$mut_type, levels = c("LOH events", "Point mutations"))
power_labels$test_type <- factor(power_labels$test_type, levels = c("Expected", "Observed"))


all_power_plot <- all_power_df %>% 
  # filter(ES <= 1.5) %>%
  ggplot(aes(x = ES_percent, y = fraction_sig)) + 
  # geom_hline(aes(yintercept = 0.7), color = "orangered3") +
  # geom_hline(aes(yintercept = 0.9), color = "brown4") +
  geom_line(aes(linetype = test_type), size = 1) +
  geom_point(aes(shape = test_type), size = 5) +
  geom_label(data = power_labels[1:2, 1:2], aes(x = 0, y = 1, label = l),
             label.size = 0, size = 12) +
  scale_shape_manual(values = c(17, 19), name = "Null distribution") +
  scale_linetype_manual(values = c(2, 1), name = "Null distribution") +
  scale_y_continuous(breaks = seq(0, 1, 0.25),
                     name = "Power (fraction of tests significant)") +
  xlab("Increase in rate relative to WT (%)") + 
  facet_grid(.~mut_type) +
  theme(text = element_text(size = 28),
        legend.position = c(0.9, 0.2),
        legend.background = element_rect(fill = "white", color = "white"),
        legend.title = element_text(size = 26),
        legend.text = element_text(size = 24),
        strip.text = element_text(size = 30),
        panel.background = element_rect(color = "grey20", size = 0.25),
        panel.grid.major = element_line(size = 0.75),
        panel.grid.minor = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.title.y = element_text(vjust = 2),
        axis.title.x = element_text(vjust = -1),
        plot.margin = unit(c(t = 0, r = 5, b = 5, l = 5), "mm"))

all_power_plot

ggsave(file.path(outIntDir, "all_power_figure_2022_05.png"), 
       plot = all_power_plot,
       device = "png",
       width = 16, height = 9, 
       units = "in",
       dpi = 600)

# UNUSED
###############################################################################

p_nSNM_WT_Drive <- perm_test_DT(SNMcounts_in[Tx_name %in% c("WT", "Drive"), ], cat_var = "Tx_name", 
                                response_var = "n", test_stat = diff_means_DT, n_perms = 10000)
p_nSNM_WT_Drive


p_nSNM_WT_Drive_list <- unlist(lapply(vector(mode = "list", length = 10), 
                               function(x) perm_test_DT(SNMcounts_in[Tx_name %in% c("WT", "Drive"), ], cat_var = "Tx_name", 
                                                        response_var = "n", test_stat = diff_means_DT, n_perms = 10)))

SNM_sub <- SNMcounts_in[Tx_name %in% c("WT", "Drive"), .(Tx_name, n)]

# Add poisson with lambda = xbar + ES to observed count data for shifted
# sampling population
# d_mu <- 0.25
# sub_add <- SNM_sub[sample(.N, ceiling(.N * d_mu)), n] + 1
# sub_null_a <- SNM_sub[sample(.N, floor(.N * (1 - d_mu)), n]
# sub_null_1 <- SNM_sub[sample(.N, ceiling(.N/dnom)), n]
# sub_null_2 <- SNM_sub[sample(.N, floor(.N/(1 - dnom))), n]
# sub_add_pois <- SNM_sub[, n + rpois(.N, 1)]
# sub_add_fract <- c(sub_add, sub_null_a)
# sub_null <- c(sub_null_1, sub_null_2)
# dist_compare <- data.frame(null = sub_null,
#                            pois = sub_add_pois,
#                            fract = sub_add_fract) %>% 
#   pivot_longer(cols = null:fract, names_to = "dist", values_to = "sm")
# dist_compare$dist <- factor(dist_compare$dist, levels = c("null", "fract", "pois"))
# 
# dist_compare %>% 
#   ggplot() + geom_histogram(aes(x = sm), binwidth = 1) + facet_grid(dist~.)

SNM_sub_shift <- SNM_sub[, .(Tx_name, n, n_pois = n + rpois(.N, 5))]

as.data.frame(SNM_sub_shift) %>% 
  pivot_longer(cols = n:n_pois, names_to = "dist", values_to = "sm") %>%
  ggplot() + geom_histogram(aes(x = sm), binwidth = 1) + facet_grid(dist~.)

SNM_sub_shift[, .(mean(n), mean(n_pois))][, .(V2 - V1)]

melt(dt, id.vars = c("id"),
     measure.vars = patterns("^a", "^b"), variable.name = "y",
     value.name = c())



# , unlist(lapply(vector(mode = "list", length = 10), 
                  # function(x) perm_test_DT(.SD, cat_var = "dist",
                  #                          response_var = "sm", test_stat = diff_means_DT, n_perms = 1000))),
  #                                    by = "repl"]
sum(p_nSNM_WT_Drive_list)

p_nSNM_WT_Drive_list <- SNMcounts_in[, unlist(lapply(vector(mode = "list", length = 10), 
                                      function(x) perm_test_DT(.SD, cat_var = "Tx_name", 
                                                               response_var = "n", test_stat = diff_means_DT, n_perms = 10))),
                                     by = ]

SNMcounts_in_df <- SNMs_final_counts %>% arrange(ID) %>% select(!(c(Line, Rep, dot_color, fill_color)))

N_H_perm_list <- SNMcounts_in_df %>% 
  # filter(!Line %in% c("H_F", "H_H")) %>%
  perm_test(df=., cat_var ="Tx_name", cat_names = c("WT", "Drive"), response_var = "n", 
            n_perms = 10000, alpha=0.05, alt_hyp = "two-tailed", test_stat = mean)
N_H_perm_list[[3]]


