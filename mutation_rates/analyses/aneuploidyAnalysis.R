# Aneuploidy analysis

qbinom(p = c(0.05, 0.95), size = 20, prob = 0.5)
qbinom(p = c(0.05, 0.95), size = 20, prob = 0.33)
qbinom(p = c(0.05, 0.95), size = 20, prob = 0.66)

pbinom(q = 3, size = 10, prob = 0.5)
pbinom(1, size = 23, prob = 0.33)

SNPs_merge_finalGT %>% 
  # filter(GT == "0/0") %>%
  # filter(Alt_DP_final >= 1) %>%
  filter(POSi %in% H_C_anc_SNPs) %>%
  filter(ID == "H_C03", CHROM == "07") %>% # pull(f_Alt) %>% min()
  ggplot() + 
  geom_line(aes(x = POSi, y = f_Alt, group = ID)) + ylim(0, 1)
  geom_histogram(aes(x = f_Alt), binwidth = 0.01) + xlim(0, 1)

nb_af <- rbinom(2831, size = 26, p = 0.346)/26


SNPs_merge_finalGT %>% summarize(s_REF = mean(Ref_DP_final), s_TOT = mean(Sum_DP_final), 
            m_f_Alt = mean(f_Alt), se_f_Alt = sd(f_Alt)/sqrt(n()),
            lo_CI = (mean(f_Alt) - sd(f_Alt)/sqrt(n())*1.96), 
            up_CI = (mean(f_Alt) + sd(f_Alt)/sqrt(n())*1.96))

dbinom(42273, size = 84546, prob = 0.66)
dbinom(138, size = 309, prob = 0.5)
dbinom(138, size = 309, prob = 0.3)

qbinom(p = c(0.01, 0.99), size = 30, prob = 0.65)

# Method for detecting aneuploidy in clones ------
# Read depth varies systematically among chromosomes and within chromosomes,
# so detecting aneuploidy by comparing mean chromosomal depth to mean genomic
# depth may be misleading. Instead:
# Divide genome into 1/3 chromosome segments
# Calculate the proportion of mean read depth at each segment (depth coefficient)
# Calculate average depth coefficient across clones at each segment. This is
# the expected relative depth at each segment.
# Calculate the ratio between the observed relative depth and expected relative depth.
# If that ratio ~1 for any 1/3 segment, no aneuploidy. 
# If it is ~0.5 for all 3 segments, this indicates the deletion of one
# homologous chromosome, if it is ~1.5 for 3/3 segments, amplification of one homolog. 

all_chr_DP <- LOH_SNPs %>% 
  # filter(GT == "0/1") %>%
  # filter(Line == "H_C") %>%
  dplyr::group_by(Tx, Line, Rep, ID, CHROM) %>% 
  dplyr::summarise(Ref_DP = mean(Ref_DP_final, na.rm = T),
                   Alt_DP = mean(Alt_DP_final, na.rm = T),
                   mean_DP = mean(Sum_DP_final, na.rm = T),
                   f_Het = sum(GT == "0/1")/sum(GT != "./."),
                   # sd_DP = sd(Sum_DP_final, na.rm = T),
                   # sd_Ref_DP = sd(Ref_DP_final, na.rm = T),
                   # f_Ref = mean(1 - f_Alt, na.rm = T),
                   # sd_f_Ref = sd(1 - f_Alt, na.rm = T),
                   # lo_CI = mean(1 - f_Alt) - sd(1 - f_Alt)/sqrt(n())*2.576, 
                   # up_CI = mean(1 - f_Alt) + sd(1 - f_Alt)/sqrt(n())*2.576,
                   N = n())

all_ID_DP <- LOH_SNPs %>% 
  # filter(GT == "0/1") %>%
  # filter(Line == "H_C") %>%
  dplyr::group_by(Tx, Line, Rep, ID) %>% 
  dplyr::summarise(Ref_DP = mean(Ref_DP_final, na.rm = T),
                   Alt_DP = mean(Alt_DP_final, na.rm = T),
                   GM_DP = mean(Sum_DP_final, na.rm = T))


# Calculate the fractional difference between the depth of each Chr in each clone
# and the mean depth of that clone
all_chr_DP <- all_chr_DP %>% 
  group_by(ID) %>% 
  mutate(f_Ref_DP = Ref_DP/mean(mean_DP),
         f_Alt_DP = Alt_DP/mean(mean_DP),
         f_Sum_DP = mean_DP/mean(mean_DP))

all_chr_DP$f_GM_DP <- operate_by_factor_match(all_ID_DP[, c("ID", "GM_DP")], 
                                              all_chr_DP[, c("ID", "mean_DP")],
                                              .fun = function(x, y) y/x)

# Calculate the fractional difference between the fractional depth of each Chr in each clone
# and the mean depth of that chromosome
all_chr_DP <- all_chr_DP %>% group_by(CHROM) %>% 
  mutate(f_CHROM_Ref_DP = f_Ref_DP/mean(f_Ref_DP),
         f_CHROM_Alt_DP = f_Alt_DP/mean(f_Alt_DP),
         f_CHROM_Sum_DP = f_Sum_DP/mean(f_Sum_DP),
         z_CHROM_Sum_DP = (f_GM_DP - mean(f_GM_DP))/sd(f_GM_DP))

z_min <- all_chr_DP %>%
  select(ID, CHROM, contains("DP")) %>%
  arrange(z_CHROM_Sum_DP) %>%
  head(n = floor(nrow(all_chr_DP) * 0.01)) %>% 
  pull(z_CHROM_Sum_DP) %>% max()

z_max <- all_chr_DP %>%
  select(ID, CHROM, contains("DP")) %>%
  arrange(desc(z_CHROM_Sum_DP)) %>%
  head(n = floor(nrow(all_chr_DP) * 0.01)) %>% 
  pull(z_CHROM_Sum_DP) %>% min()

all_chr_DP <- all_chr_DP %>% ungroup() %>% 
  mutate(f_DP = f_CHROM_Ref_DP/(f_CHROM_Ref_DP + f_CHROM_Alt_DP))

all_chr_DP %>% ggplot() + geom_histogram(aes(x = f_DP))

all_chr_DP %>% 
  # filter(f_CHROM_Alt_DP > 1.8) %>%
  filter(f_CHROM_Alt_DP > 1.8, f_CHROM_Ref_DP > 0.7) %>%
  select(ID, CHROM, f_Het, contains("DP"))

all_chr_DP %>% 
  ggplot() + 
  # geom_hline(yintercept = z_min) +
  # geom_hline(yintercept = z_max) +
  geom_line(aes(x = CHROM, y = z_CHROM_Sum_DP, group = ID))

all_chr_DP %>% 
  ggplot() + 
  geom_point(aes(x = f_CHROM_Ref_DP, y = f_CHROM_Alt_DP, color = f_Het)) +
  scale_color_viridis_b(n.breaks = 10)

all_chr_DP %>% 
  ggplot() + 
  geom_point(aes(x = z_CHROM_Sum_DP, y = f_CHROM_Ref_DP, color = f_CHROM_Alt_DP)) +
  scale_color_viridis_b(n.breaks = 10)

all_chr_DP %>% 
  ggplot() + 
  geom_point(aes(x = Ref_DP, y = Alt_DP)) +
  scale_color_viridis_b(n.breaks = 10)

sigma_r_Ref_DP <- sqrt(var(all_chr_DP$f_DP))
qnorm(0.05, mean = 0.5, sd = sigma_r_Ref_DP)

nrow(all_chr_DP)
nrow(all_chr_DP %>% filter(f_CHROM_Ref_DP > 1 - 1.96*sigma_r_Ref_DP & f_CHROM_Ref_DP < 1 + 1.96*sigma_r_Ref_DP))


all_chr_DP %>% ggplot() + geom_histogram(aes(x = f_Ref), binwidth = 0.007) + scale_y_log10()

all_chr_DP %>% 
  filter(N >= 20) %>%
  filter(f_Ref < 0.4 | f_Ref > 0.6) %>% 
  ggplot() + 
  geom_vline(aes(xintercept = 0.3333)) + 
  geom_vline(aes(xintercept = 1 - 0.3333)) + 
  geom_segment(aes(x = lo_CI, xend = up_CI, y = ID, yend = ID))

anu_DP <- SNPs_merge_finalGT %>% 
  # filter(ID == "N_F03") %>%
  # filter(GT == "0/1") %>%
  group_by(Tx_name, Line, ID, CHROM) %>% 
  summarize(mean_DP = mean(Sum_DP_final, na.rm = T),
            Ref_DP_2 = mean(Ref_DP_final, na.rm = T) * 2)

anu_DP %>% filter(Line == "N_A") %>% 
  ggplot(aes(x = CHROM)) + geom_line(aes(y = mean_DP, color = ID, group = ID))
  

# Determiniation of aneuploids:
# If there is loss of one homolog, but the other remains intact,
# we expect relative depth to be ~0 for the lost homolog,
# ~1 for the intact, and a complete LOH. If there is amplification of one homolog and
# the other remains intact, we expect ~2 and ~1 for their relative
# depths respectively, and a high proportion of heterozygosity. 
# If there are gene conversions, we expect an inverse
# relationship between relative depths with slope = -1, up to the points 
# (0,2) and (2,0), which indicate loss and amplification. If there is
# relative depth ~2 and ~0.5, it probably indicates extensive LOH and
# possible CNV. 

all_chr_DP <- LOH_SNPs %>%
  dplyr::group_by(Tx, Line, Rep, ID, CHROM) %>% 
  dplyr::summarise(Ref_DP = mean(Ref_DP_final, na.rm = T),
                   Alt_DP = mean(Alt_DP_final, na.rm = T),
                   mean_DP = mean(Sum_DP_final, na.rm = T),
                   f_Het = sum(GT == "0/1")/sum(GT != "./."))

all_chr_DP <- all_chr_DP %>% 
  group_by(ID) %>% 
  mutate(f_Ref_DP = Ref_DP/mean(mean_DP),
         f_Alt_DP = Alt_DP/mean(mean_DP),
         f_Sum_DP = mean_DP/mean(mean_DP))


# Calculate the fractional difference between the fractional depth of each Chr in each clone
# and the mean depth of that chromosome
all_chr_DP <- all_chr_DP %>% group_by(CHROM) %>% 
  mutate(f_CHROM_Ref_DP = f_Ref_DP/mean(f_Ref_DP),
         f_CHROM_Alt_DP = f_Alt_DP/mean(f_Alt_DP),
         f_CHROM_Sum_DP = f_Sum_DP/mean(f_Sum_DP))


ID_CHROM_aneuploid <- all_chr_DP %>%
  filter(((f_CHROM_Ref_DP > 1.8 | f_CHROM_Alt_DP > 1.8) & f_Het > 0.9) | 
           ((f_CHROM_Ref_DP < 0.2 | f_CHROM_Alt_DP < 0.2) & f_Het < 0.1)) %>%
  # select(ID, f_Het, contains("DP"))
  mutate(ID_CHROM = paste0(ID, "_", CHROM))


###############################################################################
all_g_DP <- SNPs_merge_finalGT %>% 
  # filter(GT == "0/1") %>%
  # filter(Line == "H_C") %>%
  dplyr::group_by(Tx, Line, Rep, ID) %>% 
  dplyr::summarise(Ref_DP = mean(Ref_DP_final, na.rm = T),
                   Alt_DP = mean(Alt_DP_final, na.rm = T),
                   mean_DP = mean(Sum_DP_final, na.rm = T),
                   sd_DP = sd(Sum_DP_final, na.rm = T),
                   sd_Ref_DP = sd(Ref_DP_final, na.rm = T),
                   f_Ref = mean(1 - f_Alt, na.rm = T),
                   sd_f_Ref = sd(1 - f_Alt, na.rm = T),
                   lo_CI = mean(1 - f_Alt) - sd(1 - f_Alt)/sqrt(n())*2.576, 
                   up_CI = mean(1 - f_Alt) + sd(1 - f_Alt)/sqrt(n())*2.576,
                   N = n())

# Windows equal to 1/3 chromosomes
subCHROMwin <- data.frame(NULL)
for (ch in chrom_bound_BY$CHROM) {
  sub_ch <- paste0(ch, ".", 1:3)
  ch_start <- chrom_bound_BY$Start[chrom_bound_BY$CHROM == ch]
  ch_end <- chrom_bound_BY$End[chrom_bound_BY$CHROM == ch]
  start_vals <- c(ch_start, (2*ch_start + ch_end)/3, (ch_start + 2*ch_end)/3)
  end_vals <- c(start_vals[c(2,3)] - 1, ch_end)
  ch_list <- list(Start = start_vals, End = end_vals, CHROM3 = sub_ch)
  subCHROMwin <- rbind(subCHROMwin, ch_list)
}

all_subCHROM_DP <- data.frame(NULL)
for (id in levels(SNPs_merge_finalGT$ID)) {
  # id = "H_C03"
  id_df <- subset(SNPs_merge_finalGT, ID == id) %>% 
    filter(GT == "0/1") %>%
    filter(POSi <= 7687467 | POSi >= 7737467)
  tx <- id_df[1, "Tx_name"]
  id_lngths <- id_df %>% 
                dplyr::group_by(CHROM) %>% 
                dplyr::summarise(chr_min = min(POS), 
                                 chr_max = max(POS))
  id_lngths$chr_ln <- id_lngths$chr_max - id_lngths$chr_min
    DP_grand_mean <- mean(id_df$Sum_DP_final, na.rm = T)
    Ref_slide_mean <- SliderCalc(df = id_df, data_col = "Ref_DP_final", index_col = "POSi",
                              window_size = subCHROMwin, summary_stat = mean)
    # Ref_slide_sd <- SliderCalc(df = id_df, data_col = "Ref_DP", index_col = "POSi",
    #                           window_size = subCHROMwin, summary_stat = sd)
    Alt_slide_mean <- SliderCalc(df = id_df, data_col = "Alt_DP_final", index_col = "POSi",
                              window_size = subCHROMwin, summary_stat = mean)
    # Alt_slide_sd <- SliderCalc(df = id_df, data_col = "Alt_DP", index_col = "POSi",
    #                           window_size = subCHROMwin, summary_stat = sd)
    sum_slide_mean <- SliderCalc(df = id_df, data_col="Sum_DP_final", index_col="POSi",
                               window_size=subCHROMwin, summary_stat = mean)
    sum_slide_sd <- SliderCalc(df = id_df, data_col="Sum_DP_final", index_col="POSi",
                               window_size=subCHROMwin, summary_stat = sd)
    slide_df <- data.frame(ID = id, subCHROM = subCHROMwin$CHROM3, 
                           DP_mean = DP_grand_mean,
                           start = Ref_slide_mean$start, end = Ref_slide_mean$end,
                           n_POSi = Ref_slide_mean$n_elements,
                           Ref_DP = Ref_slide_mean$Ref_DP_final, Alt_DP = Alt_slide_mean$Alt_DP_final, 
                           Sum_DP = sum_slide_mean$Sum_DP_final, Sum_sd = sum_slide_sd$Sum_DP_final)
    slide_df[, c("Tx_name",  "Line", "Rep")] <- c(id_df[1, c("Tx_name", "Line", "Rep")])
    all_subCHROM_DP <- rbind(all_subCHROM_DP, slide_df)
    print(id)
}
all_subCHROM_DP[, c("Tx_name",  "Line", "Rep")] <- lapply(all_subCHROM_DP[, c("Tx_name",  "Line", "Rep")], factor)


all_subCHROM_relDP <- all_subCHROM_DP %>% mutate(f_Ref = Ref_DP/Sum_DP) %>% select(-c("Ref_DP", "Alt_DP"))
all_subCHROM_relDP$CHROM <- factor(substr(all_subCHROM_relDP$subCHROM, 1, 2))


all_subCHROM_relDP$Rel_DP <- 0
all_subCHROM_relDP$ID <- factor(all_subCHROM_relDP$ID)
for (id in levels(all_subCHROM_relDP$ID)) {
  all_subCHROM_relDP$Rel_DP[all_subCHROM_relDP$ID == id] <- all_subCHROM_relDP$Sum_DP[all_subCHROM_relDP$ID == id] / 
    all_g_DP$mean_DP[all_g_DP$ID == id]
}

mean_subCHROM_relDP <- all_subCHROM_relDP %>% 
  group_by(subCHROM) %>%
  summarize(mean_Rel = mean(Rel_DP, na.rm = T),
                   sd_Rel = sd(Rel_DP, na.rm = T))

all_subCHROM_relDP$fract_Rel_DP <- 0
for (id in levels(all_subCHROM_relDP$ID)) {
  fract_rel_vals <- all_subCHROM_relDP$Rel_DP[all_subCHROM_relDP$ID == id] / mean_subCHROM_relDP$mean_Rel
  all_subCHROM_relDP$fract_Rel_DP[all_subCHROM_relDP$ID == id] <- fract_rel_vals
}

all_subCHROM_relDP$subAneu <- 2
all_subCHROM_relDP$subAneu[all_subCHROM_relDP$fract_Rel_DP == 0] <- 0
all_subCHROM_relDP$subAneu[all_subCHROM_relDP$fract_Rel_DP <= 0.6 & all_subCHROM_relDP$fract_Rel_DP > 0] <- 1
all_subCHROM_relDP$subAneu[all_subCHROM_relDP$fract_Rel_DP > 1.4 & all_subCHROM_relDP$fract_Rel_DP <= 1.9] <- 3
all_subCHROM_relDP$subAneu[all_subCHROM_relDP$fract_Rel_DP > 1.9] <- 4

all_subCHROM_relDP$aneu <- 2
for(id in levels(all_subCHROM_relDP$ID)) {
  id_idx <- all_subCHROM_relDP$ID == id
  for(ch in levels(all_subCHROM_relDP$CHROM)) {
    ch_idx <- all_subCHROM_relDP$CHROM == ch
    aneu_vals <- all_subCHROM_relDP$subAneu[id_idx & ch_idx]
    if (aneu_vals[1] == aneu_vals[2] & aneu_vals[1] == aneu_vals[3]) {
      all_subCHROM_relDP$aneu[id_idx & ch_idx] <- aneu_vals
    }
  }
}

ID_CHROM_alleleImb <- all_subCHROM_relDP %>% group_by(ID, CHROM) %>% 
  mutate(CHROM_f_Ref = mean(f_Ref), CHROM_n_POSi = sum(n_POSi)) %>% 
  distinct(ID, CHROM, .keep_all = T) %>%
  filter(CHROM_n_POSi > 30) %>% filter(CHROM_f_Ref < 0.4 | CHROM_f_Ref > 0.6) %>%
  mutate(ID_CHROM = paste0(ID, "_", CHROM)) %>% select(ID_CHROM)

ID_CHROM_aneuploid <- all_subCHROM_relDP %>% group_by(ID, CHROM) %>% 
  mutate(CHROM_n_POSi = sum(n_POSi)) %>% 
  distinct(ID, CHROM, .keep_all = T) %>%
  filter(CHROM_n_POSi > 30) %>% filter(aneu != 2) %>%
  mutate(ID_CHROM = paste0(ID, "_", CHROM)) %>% select(ID_CHROM)


all_aneuState <- all_subCHROM_relDP %>% 
                  dplyr::group_by(ID) %>%
                  distinct(CHROM, .keep_all = T) %>% 
                  select(-c(start, end, subCHROM, Rel_DP, fract_Rel_DP, subAneu))

all_aneuState$ID_CHROM <- paste0(all_aneuState$ID, "_", all_aneuState$CHROM)

aneuCounts <- all_aneuState %>%
  filter(aneu != 2) %>%
  dplyr::ungroup() %>%
  dplyr::count(Tx_name, CHROM, aneu) %>%
  complete(., Tx_name, CHROM, aneu)

aneuCounts_Tx <- all_aneuState %>%
  filter(aneu != 2) %>%
  dplyr::ungroup() %>%
  dplyr::count(Tx_name, aneu) %>%
  complete(., Tx_name, aneu)

aneuCounts$n[is.na(aneuCounts$n)] <- 0

aneuCounts$class <- factor(paste0(aneuCounts$aneu, "n"))

n_lim <- max(aneuCounts$n)

aneuCounts_wide <- aneuCounts %>% select(!c(aneu)) %>% pivot_wider(names_from = class, values_from = n)

aneu_plot <- aneuCounts %>%
  #filter(aneu != 2) %>%
  ggplot() +
  annotate(geom="rect", xmin=seq(from = 2, to = 16, by = 2) - 0.5, 
                      xmax=seq(from = 2, to = 16, by = 2) + 0.5, ymin=-Inf, ymax=Inf, alpha=0.1) + 
  geom_col(aes(x = CHROM, y = n * (as.numeric(aneu) - 2), group = Tx_name), fill = "grey70", width = 0.3,
                                 position = position_dodge(width = 1)) +
  geom_point(aes(x = CHROM, y = n * (as.numeric(aneu) - 2), group = Tx_name, color = Tx_name), size = 3,
             position = position_dodge(width = 1)) +
  # geom_segment(aes(x = CHROM, xend = CHROM, y = 0, yend = n * (as.numeric(aneu) - 2), group = Tx_name), 
  #              color = "grey50", size = 1,
  #              position = position_dodge(width = 1)) +
  geom_hline(aes(yintercept = 0), color = "white") +
  geom_hline(aes(yintercept = 0), color = "grey20", size = 0.25) +
  scale_fill_manual(values = txPal, name = "Drive Type") +
  scale_color_manual(values = txPal, name = "Drive Type") +
  scale_x_discrete(limits = rev(levels(aneuCounts$CHROM)),
                   breaks = rev(levels(aneuCounts$CHROM)),
                   labels = rev(levels(aneuCounts$CHROM)), 
                   drop = F) +
  scale_y_continuous(limits = c(-n_lim, n_lim),
                     breaks = seq(from = -n_lim, to = n_lim, by = 2),
                     labels = c(seq(from = n_lim, to = 2, by = -2), 
                                seq(from = 0, to = n_lim, by = 2))) +
  ylab("1n                            Aneuploidy Count                                3n") +
  theme(panel.grid.major.y = element_blank(), 
        panel.grid.minor.y = element_blank()) +
  coord_flip()

aneu_plot

ggsave(file.path(outIntDir, "aneu_plot_v2.png"), 
       plot = aneu_plot,
       device = "png",
       width = 11, height = 8.5, 
       units = "in",
       dpi = 300)


all_g_DP_segMeans <- all_subCHROM_DP %>% dplyr::ungroup() %>%
  dplyr::group_by(Tx_name, Line, ID) %>% 
  dplyr::summarise(mean_DP = mean(Sum_DP, na.rm = T),
                   se_DP = sd(Sum_DP, na.rm = T)/sqrt(n()))

all_aneuCount_id <- all_aneuState %>% dplyr::group_by(Tx_name, Line, ID) %>% dplyr::count(aneu) %>% filter(aneu !=2)
ids_in <- all_aneuCount_id %>% filter(Tx_name != "Drive") %>% pull(ID) %>% as.character()

subChrom_DP_in <- all_subCHROM_DP %>% dplyr::ungroup() %>%
  filter(ID %in% ids_in)

dp_lim2 <- 5*round(max(c(subChrom_DP_in$Sum_DP))/5)

# g_DP_in <- all_g_DP  %>% dplyr::ungroup() %>%
#   filter(ID %in% ids_in)

g_DP_in <- all_g_DP_segMeans %>% dplyr::ungroup() %>%
  filter(ID %in% ids_in)

aneu_call <- all_aneuState %>% 
  filter(ID %in% ids_in)

subChrom_DP_in$aneu <- factor(paste0(rep(aneu_call$aneu, each = 3), "n"))

subChrom_plot_out <- subChrom_DP_in %>% 
  ggplot() + 
  annotate(geom = "rect", xmin=chrom_bound$Start[c(TRUE, FALSE)],
           xmax=chrom_bound$End[c(TRUE, FALSE)], ymin=-Inf, ymax=Inf, alpha=0.07) +
  geom_rect(aes(xmin = start, xmax = end,
               ymin = Sum_DP - Sum_sd, ymax = Sum_DP + Sum_sd),
           fill = "grey60", alpha = 0.3) +
  geom_segment(data = g_DP_in, 
               aes(x = 0, xend = g_length,
                   y = mean_DP, yend = mean_DP),
               color = "orangered4", alpha = 0.8,
               inherit.aes = F) +
  geom_rect(data = g_DP_in,
            aes(xmin = 0, xmax = g_length,
                ymin = mean_DP - se_DP, ymax = mean_DP + se_DP),
            fill = "orangered4", alpha = 0.2, inherit.aes = F) +
  geom_label(data = g_DP_in, aes(x = 0, y = mean_DP, 
                                 label = round(mean_DP, 0)),
             label.size = 0, hjust = 1, size = 3,
             color = "orangered4") +
  geom_segment(aes(x = start, xend = end,
                   y = Sum_DP, yend = Sum_DP, color = aneu),
               size = 0.7,
               # inherit.aes = F
               ) +
  # geom_segment(aes(x = start, xend = end,
  #                  y = BY_DP, yend = BY_DP, group = label),
  #              color = "orange1",
  #              inherit.aes = F) +
  # geom_segment(aes(x = start, xend = end,
  #                  y = RM_DP, yend = RM_DP, group = label),
  #              color = "blue2",
  #              inherit.aes = F) +
  scale_x_continuous(labels=as.character(chrom_bound$CHROM),
                     breaks=(chrom_bound$Start+chrom_bound$End)/2, 
                     limits = c(0, g_length)) +
  scale_color_manual(values = allelePal) +
  xlab("Chromosome") + ylab("Mean Read Depth") +
  ylim(c(0, NA)) +
  theme(axis.text.x = element_text(size=10), 
        panel.grid.major.y = element_line(color = "grey40", size = 0.1),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.background = element_rect(fill = "white", color = "gray70"), 
        strip.background = element_rect(fill = "white", color = "gray70")) +
  facet_grid(ID~., scales = "free_y")

subChrom_plot_out


# Chi-squared test of differences -----####################
all_aneuState
aneuTable_Tx <- all_aneuState %>% ungroup() %>% mutate(class = paste0(aneu, "n")) %>% select(Tx_name, class) %>% table()

aneu_Tx_chi <- chisq.test(t(aneuTable_Tx))


M <- as.table(rbind(c(762, 327, 468), c(484, 239, 477)))
dimnames(M) <- list(gender = c("F", "M"),
                    party = c("Democrat","Independent", "Republican"))
(Xsq <- chisq.test(M))  # Prints test summary



# UNSUSED ------

# all_chrom_DPs$chrom_n <- as.numeric(all_chrom_DPs$CHROM)

chrom_DP_in <- all_chrom_DPs %>%
  # filter(ID == "N_A05")
  filter(Line == "N_A")

g_DP_in <- all_g_DP %>%
  # filter(ID == "N_A05")
  filter(Line == "N_A")


dp_lim2 <- 5*round(max(c(chrom_DP_in$BY_DP, chrom_DP_in$RM_DP))/5)

chrom_plot_out <- chrom_DP_in %>% 
  ggplot() + 
  annotate(geom="rect", xmin=chrom_bound$Start[c(TRUE, FALSE)],
           xmax=chrom_bound$End[c(TRUE, FALSE)], ymin=-Inf, ymax=Inf, alpha=0.07) +
  geom_segment(data = g_DP_in, aes(x = 1, xend = g_length,
                                   y = Sum_DP/2, yend = Sum_DP/2),
               inherit.aes = F) +
  geom_segment(aes(x = startPOSi, xend = endPOSi,
                   y = BY_DP, yend = BY_DP, group = CHROM),
               color = "orange1",
               inherit.aes = F) +
  geom_segment(aes(x = startPOSi, xend = endPOSi,
                   y = RM_DP, yend = RM_DP, group = CHROM),
               color = "blue2",
               inherit.aes = F) +
  scale_x_continuous(labels=as.character(chrom_bound$CHROM),
                     breaks=(chrom_bound$Start+chrom_bound$End)/2) +
  xlab("Chromosome") + ylab("Mean Read Depth") +
  ylim(c(0, NA)) +
  theme(axis.text.x = element_text(size=10), 
        # panel.grid.major.y = element_line(color = "grey40", size = 0.1),
        #panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.background = element_rect(fill = "white", color = "gray70"), 
        strip.background = element_rect(fill = "white", color = "gray70")) +
  facet_grid(ID~., scales = "free_y")

chrom_plot_out

# Calculates sliding window depths of allele reads and total reads for each evolved clone at sites
# initially heterozygous in the founder
all_slide_DP <- data.frame(NULL)
wndw <- 50000
for (id in levels(all_manySites$ID)) {
  # id = "N_A05"
  id_df <- subset(all_manySites, ID == id) %>% 
    # filter(GT == "0/1") %>% 
    filter(POSi <= 7687467 | POSi >= 7737467)
  tx <- id_df[1, "Tx"]
  
  id_lngths <- id_df %>% dplyr::group_by(CHROM) %>% dplyr::summarise(chr_max = max(POS), chr_min = min(POS))
  id_lngths$chr_ln <- id_lngths$chr_max - id_lngths$chr_min
  if (min(id_lngths$chr_ln) > wndw) {
    BY_slide_df <- SliderCalc(df = id_df, data_col = "BY_DP", index_col = "POS", factor_col = "CHROM", 
                              window_size = wndw, slide_interval = wndw)
    RM_slide_df <- SliderCalc(df = id_df, data_col = "RM_DP", index_col = "POS", factor_col = "CHROM", 
                              window_size = wndw, slide_interval = wndw)
    sum_slide_df <- SliderCalc(df = id_df, data_col="Sum_DP", index_col="POS", factor_col="CHROM", 
                               window_size=wndw, slide_interval = wndw)
    slide_df <- data.frame(ID = id, CHROM = BY_slide_df$CHROM, POS = BY_slide_df$start, 
                           BY_DP = BY_slide_df$value, RM_DP = RM_slide_df$value, 
                           Sum_DP = sum_slide_df$value, mean_DP = mean(sum_slide_df$value, na.rm = T))
    slide_df$CHROM <- factor(slide_df$CHROM)
    slide_df$endPOS <- slide_df$POS + wndw - 1
    slide_df$POSi <- 0
    for(ch in seq_along(levels(slide_df$CHROM))) {
      # ch = 1
      chrm_i <- slide_df$CHROM == levels(slide_df$CHROM)[ch]
      slide_df$POSi[chrm_i] <- slide_df$POS[chrm_i] + chrom_bound[ch, 1] - 1
      slide_df$endPOSi[chrm_i] <- slide_df$endPOS[chrm_i] + chrom_bound[ch, 1] - 1
    }
    slide_df$chrom_n <- as.numeric(slide_df$CHROM)
    slide_df[, c("Tx",  "Line", "Rep")] <- c(id_df[1, c("Tx", "Line", "Rep")])
    all_slide_DP <- rbind(all_slide_DP, slide_df)
    print(id)
  }
  all_slide_DP[, c("Tx",  "Line", "Rep")] <- lapply(all_slide_DP[, c("Tx",  "Line", "Rep")], factor)
}



slide_plot_in <- all_slide_DP %>% 
  filter(Line == "H_A")
# filter(Tx == "H")

id_meanDP <- slide_plot_in %>% 
  dplyr::ungroup() %>%
  dplyr::group_by(Rep) %>% 
  dplyr::summarize(minPOSi = min(POSi),
                   maxPOSi = max(POSi), 
                   BY_DP = mean(BY_DP, na.rm = T), 
                   RM_DP = mean(RM_DP, na.rm = T), 
                   Sum_DP = mean(Sum_DP, na.rm = T))

dp_lim1 <- 5*ceiling(max(c(slide_plot_in$BY_DP, slide_plot_in$RM_DP), na.rm = T)/5)

slide_plot_out <- slide_plot_in %>% 
  ggplot() + 
  annotate(geom="rect", xmin=chrom_bound$Start[c(TRUE, FALSE)],
           xmax=chrom_bound$End[c(TRUE, FALSE)], ymin=-Inf, ymax=Inf, alpha=0.07) +
  geom_line(aes(x = POSi, y = BY_DP, group = CHROM), color = "orange1", size = 0.5) +
  geom_line(aes(x = POSi, y = RM_DP, group = CHROM), color = "blue2", size = 0.5) +
  geom_line(aes(x = POSi, y = Sum_DP, group = CHROM), color = "grey50", size = 0.5) +
  # geom_segment(data = id_meanDP, aes(x = minPOSi, xend = maxPOSi, y = Sum_DP, yend = Sum_DP), color = "grey50", size = 0.5) +
  scale_x_continuous(labels=as.character(chrom_bound$CHROM),
                     breaks=(chrom_bound$Start+chrom_bound$End)/2) +
  xlab("Chromosome") + ylab("Mean Read Depth") +
  # ylim(c(0, dp_lim1)) +
  theme(axis.text.x = element_text(size=10), 
        # panel.grid.major.y = element_line(color = "grey40", size = 0.1),
        #panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.background = element_rect(fill = "white", color = "gray70"), 
        strip.background = element_rect(fill = "white", color = "gray70")) +
  facet_grid(Rep~., scales = "free_y")

slide_plot_out

# Chromosome mean depths to detect aneuploidy events

all_chrom_DPs <- all_manySites %>% 
  dplyr::group_by(Tx, Line, Rep, ID, CHROM) %>% 
  dplyr::summarise(BY_DP = mean(BY_DP, na.rm = T),
                   RM_DP = mean(RM_DP, na.rm = T),
                   Sum_DP = mean(Sum_DP, na.rm = T))

all_chrom_DPs$startPOSi <- 0
all_chrom_DPs$endPOSi <- 0
for(ch in levels(all_chrom_DPs$CHROM)) {
  all_chrom_DPs$startPOSi[all_chrom_DPs$CHROM == ch] <- chrom_bound$Start[chrom_bound$CHROM == ch]
  all_chrom_DPs$endPOSi[all_chrom_DPs$CHROM == ch] <- chrom_bound$End[chrom_bound$CHROM == ch]
}

