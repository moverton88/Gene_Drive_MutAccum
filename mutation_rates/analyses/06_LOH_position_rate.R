###############################################################################
# Number of times a SNP has been converted ####
LOH_SNPs$Tx_ID <- Recode_Tx_ID(LOH_SNPs$Tx)

n_clones_LOH_Line <- LOH_SNPs %>% filter(Rep != "00") %>% distinct(Line, ID) %>% count(Line)

SNPs_counts_Line_GT <- LOH_SNPs %>%
  filter(Rep != "00") %>%
  count(Line, CHROM, POS, GT, name = "n_conv") %>%
  ungroup %>%
  mutate(per_clone = operate_by_factor_match(n_clones_LOH_Line[, c("Line", "n")], 
                                             data.frame(Line, n_conv),
                                             function(x, y) y/x)) 
mean_m_cover_Line <- SNPs_counts_Line_GT %>% 
  group_by(Line, CHROM, POS) %>%
  summarize(n_m = sum(per_clone)) %>%
  group_by(Line) %>%
  summarize(mean_m = mean(n_m))

mean_m_cover_Line <- SNPs_counts_Line_GT %>% 
  group_by(Line, CHROM, POS) %>%
  summarize(n_m = sum(n_conv)) %>%
  ungroup() %>%
  # filter(n_m > 2) %>%
  count(Line)

SNPs_counts_Line <- LOH_SNPs %>%
  filter(Rep != "00") %>%
  count(Line, CHROM, POS, name = "n_clones") %>%
  ungroup %>%
  mutate(per_clone = operate_by_factor_match(n_clones_LOH_Line[, c("Line", "n")], 
                                             data.frame(Line, n_clones),
                                             function(x, y) y/x)) 


# n_LOH against marker density by Line
mean_m_cover_Line <- SNPs_counts_Line %>% 
  group_by(Line) %>%
  summarize(mean_m = mean(per_clone),
            mean_hi = sum(per_clone > 0.5)/n(),
            sum_m = n())

LOH_Line <- all_LOHcounts_merge_NS %>% group_by(Line) %>% summarize(tot_LOH = sum(n_LOH))

mean_m_cover_Line <- merge(mean_m_cover_Line, LOH_Line, by = "Line")

mean_m_cover_Line %>% 
  ggplot() + geom_point(aes(x = mean_m, y = tot_LOH)) +
  xlim(0, NA) +
  ylim(0, NA)

###############################################################################
# Supp Fig. Number of markers vs number of LOH by clone ####
m_cover_ID <- LOH_SNPs %>% count(ID, name = "n_m")
m_cover_ID <- merge(m_cover_ID, all_LOHcounts_merge_NS[, c("ID", "n_LOH")], by = "ID")


marker_cover_LOH_lm <- lm(n_LOH ~ n_m, data = m_cover_ID)
summary(marker_cover_LOH_lm)

mcL_coeff <- coefficients(marker_cover_LOH_lm)
marker_cover_LOH_seg <- data.frame(x.min = min(m_cover_ID$n_m),
                                   x.max = max(m_cover_ID$n_m),
                                   y.min = min(m_cover_ID$n_m) * mcL_coeff[2] + mcL_coeff[1],
                                   y.max = max(m_cover_ID$n_m) * mcL_coeff[2] + mcL_coeff[1])

m_cover_ID %>% 
  ggplot() + 
  geom_point(aes(x = n_m, y = n_LOH)) +
  geom_segment(data = marker_cover_LOH_seg, 
               aes(x = x.min, xend = x.max, y = y.min, yend = y.max)) +
  # xlim(0, NA) +
  ylim(0, NA)

###############################################################################
# Supp Fig. Marker density vs LOH rate by chromosome ####
m_cover_CHROM <- LOH_SNPs %>% 
  count(CHROM, name = "n_m") %>%
  mutate(m_rate = operate_by_factor_match(chrom_lengths_BY_df[, c("CHROM", "chrom_length")],
                                          data.frame(CHROM, n_m),
                                          function(x, y) y/x))
LOH_rate_CHROM <- all_LOHbounds_merge_NS %>% 
  count(CHROM, name = "n_LOH") %>%
  mutate(LOH_rate = operate_by_factor_match(chrom_lengths_BY_df[, c("CHROM", "chrom_length")],
                                          data.frame(CHROM, n_LOH),
                                          function(x, y) y/x))

m_cover_CHROM <- merge(m_cover_CHROM, LOH_rate_CHROM, by = "CHROM")


marker_cover_CHROM_lm <- lm(LOH_rate ~ m_rate, data = m_cover_CHROM)
summary(marker_cover_CHROM_lm)

mcL_coeff <- coefficients(marker_cover_CHROM_lm)
marker_cover_CHROM_seg <- data.frame(x.min = min(m_cover_CHROM$m_rate),
                                   x.max = max(m_cover_CHROM$m_rate),
                                   y.min = min(m_cover_CHROM$m_rate) * mcL_coeff[2] + mcL_coeff[1],
                                   y.max = max(m_cover_CHROM$m_rate) * mcL_coeff[2] + mcL_coeff[1])

m_cover_CHROM %>% 
  ggplot() + 
  geom_point(aes(x = m_rate, y = LOH_rate)) +
  geom_segment(data = marker_cover_CHROM_seg, 
               aes(x = x.min, xend = x.max, y = y.min, yend = y.max)) +
  # xlim(0, NA) +
  ylim(0, NA)

###############################################################################

# Distribution of fraction of total marker coverage 
SNPs_counts_Line %>% 
  group_by(Line, CHROM, POS) %>%
  summarize(n_m = sum(per_clone)) %>%
  ungroup() %>%
  ggplot() + 
  geom_histogram(aes(x = n_m), binwidth = 0.1) +
  xlim(-0.05, 1.05) +
  facet_wrap(~Line)

# Fraction of sites have maerker coverage > 0.25
SNPs_counts_Line %>% 
  group_by(Line, CHROM, POS) %>%
  summarize(n_m = sum(per_clone), 
            f_m = sum(per_clone) >= 0.5) %>%
  group_by(Line) %>%
  summarize(f_n_m = mean(n_m), f_f_m = sum(f_m)/n())

SNPs_counts_Line %>% 
  ggplot() + geom

SNPs_counts_Line %>% 
  select(!n_conv) %>%
  pivot_wider(names_from = "GT", values_from = "per_clone", values_fill = 0) %>%
  rename(Ref = `0/0`, het = `0/1`, Alt = `1/1`) %>%
  # filter(het >= 0.25) %>%
  mutate(LOH = Ref + Alt) %>%
  ggplot() + 
  # geom_histogram(aes(x = LOH))
  geom_jitter(aes(x = (LOH + het), y = LOH, #  + 0.1*(2 - as.numeric(Tx_ID)
                  color = Line),
              height = 0.01, width = 0.01,
              size = 0.1) +
  # geom_point(aes(x = POS, y = (LOH + het), #  + 0.1*(2 - as.numeric(Tx_ID)
  #                color = Tx_ID), size = 0.3) +
  # geom_point(aes(x = POS, y = LOH + 0.1*(2 - as.numeric(Tx_ID)), 
  #                color = Tx_ID), size = 0.5) +
  # geom_line(aes(x = POS, y = n_conv, color = Tx_ID)) +
  # scale_y_continuous(limits = c(0, 1)) +
  facet_grid(.~Line) +
  # facet_grid(Tx_ID~CHROM, scales = "free_x") +
  scale_color_manual(values = txPal)

SNPs_counts_Line %>%
  mutate(LOH = GT != "0/1") %>%
  # filter(CHROM == "05") %>%
  ggplot() + 
  geom_histogram(aes(x = n_conv), binwidth = 1) + 
  facet_grid(Line~LOH, scales = "free_y")

SNPs_counts_Tx <- LOH_SNPs %>%
  # filter(Tx_ID == "D") %>%
  # filter(CHROM == "05") %>% 
  filter(Rep != "00") %>%
  count(Tx_ID, CHROM, POS, POSi, GT, name = "n_conv") %>%
  ungroup %>%
  # group_by(Tx_ID) %>%
  # group_modify(.$n_conv, .fun = function(x, y) x/y, y = n_clones_LOH_xTx$n)
  mutate(per_clone = operate_by_factor_match(n_clones_LOH_xTx[, c("Tx_ID", "n")], 
                                             data.frame(Tx_ID, n_conv),
                                             function(x, y) y/x)) 

SNPs_counts <- LOH_SNPs %>%
  # filter(Tx_ID == "D") %>%
  # filter(CHROM == "05") %>% 
  filter(Rep != "00") %>%
  count(CHROM, POS, POSi, name = "n_clones") %>%
  ungroup %>%
  # group_by(Tx_ID) %>%
  # group_modify(.$n_conv, .fun = function(x, y) x/y, y = n_clones_LOH_xTx$n)
  mutate(per_clone = n_clones/sum(n_clones_LOH_xTx$n)) 

mean(SNPs_counts$per_clone)
nrow(SNPs_counts %>% filter(n_clones >= sum(n_clones_LOH_xTx$n)*0.75))/nrow(SNPs_counts)

# Distribution of fraction of total marker coverage 
SNPs_counts_Tx %>% 
  group_by(Tx_ID, POSi) %>%
  summarize(n_m = sum(per_clone)) %>%
  ggplot() + 
  geom_histogram(aes(x = n_m), binwidth = 0.05) + 
  facet_grid(Tx_ID~.)

SNPs_counts_Line %>%
  ggplot() + 
  geom_histogram(aes(x = per_clone), binwidth = 0.05)
  
# Fraction of sites have maerker coverage > 0.25
SNPs_counts_Tx %>% 
  group_by(POSi) %>%
  summarize(n_m = sum(per_clone), 
            f_m = sum(per_clone) > 0.5) %>%
  group_by(Tx_ID) %>%
  summarize(f_f_m = sum(f_m)/n())

SNPs_counts_Tx %>% 
  select(!n_conv) %>%
  pivot_wider(names_from = "GT", values_from = "per_clone", values_fill = 0) %>%
  rename(Ref = `0/0`, het = `0/1`, Alt = `1/1`) %>%
  # filter(het >= 0.25) %>%
  mutate(LOH = Ref + Alt) %>%
  # filter(CHROM == "09") %>%
  # filter(POS > 434000) %>%
  # mutate(LOH = GT != "0/1") %>%
  # filter(GT != "0/1") %>%
  # filter(Tx_ID == "W") %>%
  ggplot() + 
  # geom_histogram(aes(x = LOH))
  # geom_jitter(aes(x = (LOH + het), y = LOH, #  + 0.1*(2 - as.numeric(Tx_ID)
  #                # color = Tx_ID
  #                ),
  #             height = 0.01, width = 0.01,
  #            size = 0.1) +
  geom_point(aes(x = POS, y = (LOH + het), #  + 0.1*(2 - as.numeric(Tx_ID)
                 color = Tx_ID
                 ), size = 0.3) +
  # geom_point(aes(x = POS, y = LOH + 0.1*(2 - as.numeric(Tx_ID)), 
  #                color = Tx_ID), size = 0.5) +
  # geom_line(aes(x = POS, y = n_conv, color = Tx_ID)) +
  # scale_y_continuous(limits = c(0, 1)) +
  facet_grid(Tx_ID~CHROM) +
  # facet_grid(Tx_ID~CHROM, scales = "free_x") +
  scale_color_manual(values = txPal)

SNPs_counts_Tx %>%
  mutate(LOH = GT != "0/1") %>%
  # filter(CHROM == "05") %>%
  ggplot() + 
  geom_histogram(aes(x = n_conv), binwidth = 1) + 
  facet_grid(Tx_ID~LOH, scales = "free_y")

LOH_SNPs %>% 
  count(Line, ID) %>% 
  separate(ID, sep = -2, into = c("Line", "Rep"), remove = F) %>%
  ggplot() + 
  geom_boxplot(aes(x = Line, y = n)) +
  geom_point(aes(x = Line, y = n, color = Rep == "00"), alpha = 0.7) +
  scale_color_manual(values = c("grey60", "Red3")) +
  ylim(0, NA)

marker_SW <- LOH_SNPs %>%
  mutate(m = 1) %>%
  SliderCalc(data_col = "m", index_col = "POSi", 
             window_size = 50000, slide_interval = 50000,
             summary_stat = sum, factor_col = "ID", chrom_win = T)

# marker_SW_o <- marker_SW

marker_SW <- marker_SW %>% mutate(start_POSi = round(start), end_POSi = round(end)) %>% select(!start:end)
s_POS <- marker_SW %>% ConvertPosIndicies(pos_col = "start_POSi", index_out = "POS", add_chroms = T)
names(s_POS)[1] <- "start_POS"
e_POS <- marker_SW %>% ConvertPosIndicies(pos_col = "end_POSi", index_out = "POS")
marker_SW <- data.frame(marker_SW, s_POS[c(2,1)], end_POS = e_POS)
# marker_SW$Tx_ID <- factor(marker_SW$Tx_ID, levels = c("W", "C", "D"))
marker_SW <- marker_SW %>% 
  mutate(mid_POS = start_POS + 50000/2, mid_POSi = (start_POSi + end_POSi)/2,
         mid_POS_kb = (start_POS + 50000/2)/1000, mid_POSi_kb = (start_POSi + end_POSi)/2000)

marker_SW$rom_CHROM <- factor(chrom_IDs$rom_CHROM[as.numeric(marker_SW$CHROM)], levels = roman_chr)

marker_SW_plot <- marker_SW %>% 
  # filter(Line != "F_D") %>%
  # filter(CHROM == "03") %>%
  # filter(Tx_ID == "D") %>%
  ggplot(aes(x = mid_POS, y = log10(m))) + 
  geom_vline(data = centrom_df, aes(xintercept = POS/1000), 
             linetype = 2,color = "blue4", size = 0.25) +
  geom_line(aes(color = Line), show.legend = F) +
  # xlim(0, max(chrom_lengths_BY)) +
  xlab("Window Position (kb)") + ylab("Breakpoints in window / Breakpoints in chromosome") +
  xlim(0, NA) + ylim(0, NA) +
  # scale_color_manual(values = txPal, name = "Strain") + 
  # scale_fill_manual(values = txPal, name = "Strain") + 
  facet_wrap(~rom_CHROM, ncol = 4, scales = "free_x") +
  theme(
    # legend.text = element_text(size = 13),
        # legend.title = element_text(size = 14),
        strip.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 16),
        # text = element_text(size = 18), 
        plot.margin = margin(20, 20, 20, 20),
        panel.grid.minor = element_blank(),
        # legend.position = c(0.97, 0.94),
        # legend.background = element_rect(fill = "white", color = "white")
        )

marker_SW_plot

marker_SW <- merge(marker_SW, 
                   all_LOHrates_SW50k[, c("ID", "CHROM", "start_POSi", "LOH_rate")],
                   by = c("ID", "CHROM", "start_POSi"))

LOH_ID_SW <- LOH_ID_SW %>% mutate(start_POSi = start)
marker_SW <- merge(marker_SW, 
                   LOH_ID_SW[, c("ID", "CHROM", "start_POSi", "bp")],
                   by = c("ID", "CHROM", "start_POSi"))

marker_lm <- lm(LOH_rate ~ log10(m), data = marker_SW %>% filter(m > 0))
summary(marker_lm)
marker_bp_lm <- lm(bp ~ m, data = marker_SW %>% filter(m > 0))
summary(marker_bp_lm)

###############################################################################
# Supp Fig. The dependence of LOH detection on number of markers sampled in a 50kb window ####

marker_SW_means <- marker_SW %>% 
  group_by(start_POSi) %>%
  summarize(mean_m = mean(m, na.rm = T), 
            log_mean_m = log10(mean(m, na.rm = T)),
            mean_LOH = mean(LOH_rate, na.rm = T))

marker_lm <- lm(mean_LOH ~ log_mean_m, data = marker_SW_means)
summary(marker_lm)

mc_coeff <- coefficients(marker_lm)
marker_cover_seg <- data.frame(x.min = min(marker_SW_means$log_mean_m),
                                   x.max = max(marker_SW_means$log_mean_m),
                                   y.min = min(marker_SW_means$log_mean_m) * mc_coeff[2] + mc_coeff[1],
                                   y.max = max(marker_SW_means$log_mean_m) * mc_coeff[2] + mc_coeff[1])


marker_SW_means %>% ggplot() +
  geom_point(aes(x = log_mean_m, y = mean_LOH)) +
  geom_segment(data = marker_cover_seg, aes(x = x.min, xend = x.max, y = y.min, yend = y.max))

###############################################################################

marker_SW %>% ggplot() +
  geom_jitter(aes(x = m, y = bp), 
              height = 0.03, width = 0.03,
              size = 0.25)


bin_values <- function(x, binwidth = NULL, n_bins = NULL, breaks = NULL,
                       output = c("values", "counts")) {
  n_na <- sum(is.na(x)) 
  x_clean <- x[!is.na(x)]
  min_x <- min(x_clean)
  max_x <- max(x_clean)
  if(sum(c(length(binwidth) > 0, length(n_bins) > 0, length(breaks) > 0)) == 1) {
    if(length(binwidth) > 0) {
    # binwidth
      bins <- seq(min_x, max_x + binwidth, binwidth)
    } else if(length(n_bins) > 0) {
      # n_bins 
      bins <- seq(min_x, max_x, length.out = n_bins)
    } else if(length(breaks) > 0) {
      # n_bins 
      bins <- breaks
      bins <- sort(bins)
      bins <- c(bins, ifelse(min(bins) > min_x, min_x, min(bins)))
      bins <- c(bins, ifelse(max(bins) < max_x, max_x, max(bins)))
      bins <- sort(unique(bins))
      
    } else {
      print("Must provide one of binwidth, n_bins, or breaks")
    }
  } else {
    print("Must provide only one of binwidth, n_bins, or breaks")
  }
  bins_mid <- (lag(bins)[-1] + bins[-1])/2
  bins_list <- as.list(bins)
  hist_low <- lapply(bins_list, function(y) x_clean < y)
  hist_low <- hist_low[2:length(hist_low)]
  hist_high <- lapply(bins_list, function(y) x_clean >= y)
  hist_high <- hist_high[1:(length(hist_high) - 1)]
  hist_list <- list(NULL)
  for(l in 1:length(hist_low)) {
    hist_bt <- mapply(function(w, z) w & z, hist_low[l], hist_high[l])
    hist_list[[l]] <- hist_bt
  }
  hist_names <- paste0("bin_", lag(bins)[-1], "_", bins[-1])
  bins_low <- lag(bins)[-1]
  names(hist_list) <- hist_names
  hist_df <- as.data.frame(hist_list)
  hist_counts <- apply(hist_df, MARGIN = 2, function(y) sum(y))
  i_bins <- apply(hist_df, MARGIN = 1, function(y) which(y))
  names_col <- hist_names[i_bins]
  bins_col <- bins_mid[i_bins]
  hist_values <- data.frame(value = x, bin_name = names_col, bin_mid = bins_col)
  if(output == "counts") {
    return(hist_counts)
  } else if(output == "values") {
    return(hist_values)
  } else if(sum(output %in% c("values", "counts")) == 2) {
    hist_list_out <- list(counts = hist_counts, values = hist_values)
    return(hist_list_out)
  }
  
}

marker_SW$m_bin_name <- bin_values(marker_SW$m, binwidth = 10, output = "values")[, 2]
marker_SW$m_bin <- bin_values(marker_SW$m, binwidth = 10, output = "values")[, 3]

# marker_SW$m_bin <- factor(marker_SW$m_bin, levels = sort(unique(marker_SW$m_bin)))
marker_SW$LOH_bin <- ifelse(marker_SW$LOH_rate > 0, 1, 0)

marker_bin_mLOH <- marker_SW %>% group_by(m_bin) %>% summarise(m_LOH = mean(LOH_rate, na.rm = T))
marker_bin_mLOH %>% ggplot() + geom_point(aes(x = m_bin, y = m_LOH))

marker_LOH_bin <- table(marker_SW$m_bin, marker_SW$LOH_bin) %>% as.data.frame() %>% arrange(Var1)

marker_LOH_bin %>% ggplot() + geom_col(aes(x = Var1, y = Freq, fill = Var2), position = position_dodge())

marker_SW %>% ggplot() +
  geom_histogram(aes(x = log10(m)), 
              binwidth = 0.1) +                                                                    
  facet_grid(LOH_rate >= 0.5~., scales = "free_y")

marker_SW %>% ggplot() +
  geom_histogram(aes(x = LOH_rate), 
                 binwidth = 0.1) +                                                                    
  facet_grid(log10(m) >= 1~., scales = "free_y")

marker_POSi_SW <- marker_SW %>%
  filter(m > 10) %>%
  count(start_POSi, CHROM, rom_CHROM, start_POS)
  
marker_ID_SW_plot <- marker_POSi_SW %>% 
  ggplot(aes(x = start_POS, y = n)) + 
  geom_vline(data = centrom_df, aes(xintercept = POS/1000), 
             linetype = 2,color = "blue4", size = 0.25) +
  geom_line() +
  # xlim(0, max(chrom_lengths_BY)) +
  xlab("Window Position (kb)") + ylab("Number of clones with marker") +
  xlim(0, NA) + ylim(0, NA) +
  # scale_color_manual(values = txPal, name = "Strain") + 
  # scale_fill_manual(values = txPal, name = "Strain") + 
  facet_wrap(~rom_CHROM, ncol = 4, scales = "free_x") +
  theme(
    # legend.text = element_text(size = 13),
    # legend.title = element_text(size = 14),
    strip.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 16),
    # text = element_text(size = 18), 
    plot.margin = margin(20, 20, 20, 20),
    panel.grid.minor = element_blank(),
    # legend.position = c(0.97, 0.94),
    # legend.background = element_rect(fill = "white", color = "white")
  )
marker_ID_SW_plot

marker_SW_Tx <- LOH_SNPs %>%
  mutate(m = 1) %>%
  SliderCalc(data_col = "m", index_col = "POSi", 
           window_size = 50000, slide_interval = 50000,
           summary_stat = sum, factor_col = "Tx_ID", chrom_win = T)

marker_SW_Tx <- merge(marker_SW_Tx, 
                   LOH_mid_SW[, c("Tx_ID", "start", "CHROM", "bp")], 
                   by = c("Tx_ID", "start", "CHROM"))

marker_SW_Tx <- marker_SW_Tx %>% mutate(start_POSi = round(start), end_POSi = round(end)) %>% select(!start:end)
s_POS <- marker_SW_Tx %>% ConvertPosIndicies(pos_col = "start_POSi", index_out = "POS", add_chroms = T)
names(s_POS)[1] <- "start_POS"
e_POS <- marker_SW_Tx %>% ConvertPosIndicies(pos_col = "end_POSi", index_out = "POS")
marker_SW_Tx <- data.frame(marker_SW_Tx, s_POS[c(2,1)], end_POS = e_POS)
marker_SW_Tx$Tx_ID <- factor(marker_SW_Tx$Tx_ID, levels = c("W", "C", "D"))
marker_SW_Tx <- marker_SW_Tx %>% 
  mutate(mid_POS = start_POS + 50000/2, mid_POSi = (start_POSi + end_POSi)/2,
         mid_POS_kb = (start_POS + 50000/2)/1000, mid_POSi_kb = (start_POSi + end_POSi)/2000)

marker_SW_Tx$rom_CHROM <- factor(chrom_IDs$rom_CHROM[as.numeric(marker_SW_Tx$CHROM)], levels = roman_chr)

marker_SW_Tx_plot <- marker_SW_Tx %>% 
  # filter(CHROM == "03") %>%
  # filter(Tx_name == "Drive") %>%
  ggplot(aes(x = mid_POS, y = log10(m))) + 
  geom_vline(data = centrom_df, aes(xintercept = POS/1000), 
             linetype = 2,color = "blue4", size = 0.25) +
  # geom_ribbon(aes(ymin = CI_lo, ymax = CI_up, fill = Tx_name, color = Tx_name), 
  #             size = 0.05, alpha = 0.1) +
  geom_line(aes(color = Tx_ID)) +
  # geom_line(data = genome_scores_50kb, 
  #           aes(x = mid_POS/1000, y = -(score - min_score)/((max(score) - min_score)*3))) +
  # geom_line(aes(x = mid_POS_kb, y = log(rate, 10), color = Tx_name)) +
  # geom_point(aes(x = mid_POS_kb, y = rate, color = Tx_name)) +
  # xlim(0, max(chrom_lengths_BY)) +
  xlab("Window Position (kb)") + ylab("Breakpoints in window / Breakpoints in chromosome") +
  xlim(0, NA) + ylim(0, NA) +
  # scale_y_continuous(limits = c(-1/3, 0.4), 
  #                    breaks = c(-1/3, -2/9, -1/9, seq(0, 0.4, 0.1)), 
  #                    labels = c(210, 140, 70, seq(0, 0.4, 0.1))) +
  scale_color_manual(values = txPal, name = "Strain") + 
  scale_fill_manual(values = txPal, name = "Strain") + 
  facet_wrap(~rom_CHROM, ncol = 4, scales = "free_x") +
  theme(legend.text = element_text(size = 13),
        legend.title = element_text(size = 14),
        strip.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 16),
        # text = element_text(size = 18), 
        plot.margin = margin(20, 20, 20, 20),
        panel.grid.minor = element_blank(),
        legend.position = c(0.97, 0.94),
        legend.background = element_rect(fill = "white", color = "white"))

marker_SW_Tx_plot

marker_SW_Tx %>%
  ggplot() +
  geom_jitter(aes(x = m, y = bp, color = Tx_ID), height = 0.3, width = 0, size = 0.25)

marker_SW_Tx %>%
  ggplot() +
  geom_histogram(aes(x = m, fill = Tx_ID), binwidth = 1000) +
  facet_grid(Tx_ID~.)


###############################################################################
# Prepare LOH position data ------
# Calculated as n_LOH_events/n_clones/gens
POSi_data_in <- all_LOHbounds_merge_NS
  # filter(!isTerm))

POSi_het <- all_GT_bounds_merge %>% 
  filter(GT == "0/1")


POSi_data_in <- all_LOHbounds_merge_NS %>% get_LOH_coordinates()

POSi_het <- all_GT_bounds_merge %>% 
  filter(GT == "0/1") %>% get_LOH_coordinates(Het = T)

POSi_data_wHet <- merge(POSi_data_in, POSi_het, all = T) %>% arrange(ID, est_start)
POSi_data_wHet$GT <- factor(POSi_data_wHet$GT, levels = GT_levels)

LOH_bounds_plot <- POSi_data_wHet %>%
  # filter(CHROM == "09") %>%
  # filter(Tx_ID == "W") %>%
  # filter(est_end > 580000) %>%
  # filter(GT != "0/1", .preserve = T) %>% 
  # arrange(start_POSi) %>% 
  # count(length > 2)
  ggplot() + 
  # geom_segment(aes(x = est_start, xend = ifelse(est_length < 10000, est_start + 10000, est_end), 
                   # y = ID, yend = ID, color = GT), size = 0.75) +
  geom_point(aes(x = est_mid, y = ID, color = GT), size = 0.5) +
  scale_color_manual(values = c(allelePal[1], "grey80", allelePal[3]), drop = F) +
  # scale_color_manual(values = allelePal[c(1,3)]) +
  scale_x_continuous(expand = c(0,0)) +
  facet_wrap(~CHROM, scales = "free_x")

LOH_bounds_plot

ggsave(file.path(outIntDir, "LOH_bounds_plot.png"), 
       plot = LOH_bounds_plot,
       device = "png",
       width = 50, height = 30, 
       units = "in",
       dpi = 600, limitsize = FALSE)

# POSi_data_wHet %>% filter(CHROM == "09", est_end > 5792574, GT != "0/1") %>% count(Tx_ID)

# Use the LOHrateSW() function to determine clone/chromosomes too few markers
all_LOHrates_Chrm <- POSi_data_wHet %>% 
  LOHrateSW(by_GT = F, win = "CHROM", 
            chrom_df = chrom_bound_BY,
            rm_zero_valids = F)

low_cover_ID_CHROM <- all_LOHrates_Chrm %>% 
  filter(prop_valid < 0.8) %>% select(ID, CHROM) %>% 
  mutate(ID_CHROM = paste0(ID, "_", CHROM))

POSi_data_in$ID_CHROM <- paste0(POSi_data_in$ID, "_", POSi_data_in$CHROM)

POSi_data_in$lc_CHROM <- POSi_data_in$ID_CHROM %in% low_cover_ID_CHROM$ID_CHROM

###############################################################################
# Calculate LOH event rates against chromosome length #########################

LOHcounts_ID_CHROM <- POSi_data_in %>% 
  group_by(ID, CHROM, .drop = F) %>% 
    count(name = "n_LOH") %>% as.data.frame()
LOHcounts_ID_CHROM$CHROM <- factor(LOHcounts_ID_CHROM$CHROM)
LOHcounts_ID_CHROM$chrom_length <- 0
LOHcounts_ID_CHROM$chrom_length <- operate_by_factor_match(chrom_lengths_BY_df[, c("CHROM", "chrom_length")],
                          LOHcounts_ID_CHROM[, c("CHROM", "chrom_length")],
                          .fun = function(x, y) x+y)
LOHcounts_ID_CHROM <- CategoriesFromID(LOHcounts_ID_CHROM)

###############################################################################
# LOH rate for each strain across chromosomes ######
LOHcounts_ID_CHROM_mean <- LOHcounts_ID_CHROM %>% 
  group_by(Tx_ID, CHROM) %>% 
  summarize(total_LOH = sum(n_LOH), mean_LOH = mean(n_LOH), se_LOH = se(n_LOH), 
            mean_rate = sum(n_LOH)/n(), se_rate = se(n_LOH)/n(), 
            CI_rate = 1.96*se(n_LOH)/n())

LOH_CHROM_CI_list <- LOHcounts_ID_CHROM %>% 
  group_by(Tx_ID, CHROM) %>%
  group_modify(~ data.frame(boot.ci(boot(.x$n_LOH,
                                         statistic = sum.fun, R = 1000), 
                                    type = "perc")$percent[4:5] %>% t())) %>%
  rename(CI_95lo = X1, CI_95up = X2) %>% 
  select(CI_95lo, CI_95up)

LOHcounts_ID_CHROM_mean <- merge(LOHcounts_ID_CHROM_mean, 
                                   LOH_CHROM_CI_list,
                                   by = c("Tx_ID", "CHROM"))

LOHcounts_ID_CHROM_mean$chrom_length <- rep(chrom_lengths_BY, 3)
LOHcounts_ID_CHROM_mean <- LOHcounts_ID_CHROM_mean %>% 
  group_by(Tx_ID, CHROM) %>%
  mutate(bp_rate = total_LOH/chrom_length,
         bp_CI_lo = CI_95lo/chrom_length,
         bp_CI_up = CI_95up/chrom_length)

LOHcounts_ID_CHROM_mean$bp_rate <- operate_by_factor_match(
  n_clones_LOH_xTx,
  LOHcounts_ID_CHROM_mean[, c("Tx_ID", "bp_rate")],
  .fun = function(x, y) y/x/n_gens)
LOHcounts_ID_CHROM_mean$bp_CI_lo <- operate_by_factor_match(
  n_clones_LOH_xTx,
  LOHcounts_ID_CHROM_mean[, c("Tx_ID", "bp_CI_lo")],
  .fun = function(x, y) y/x/n_gens)
LOHcounts_ID_CHROM_mean$bp_CI_up <- operate_by_factor_match(
  n_clones_LOH_xTx,
  LOHcounts_ID_CHROM_mean[, c("Tx_ID", "bp_CI_up")],
  .fun = function(x, y) y/x/n_gens)

LOHcounts_ID_CHROM_mean$rom_CHROM <- 
  roman_chr[as.numeric(LOHcounts_ID_CHROM_mean$CHROM)]

aov_LOH_Tx_CHROM <- aov(bp_rate ~ CHROM + Tx_ID, data = LOHcounts_ID_CHROM)
summary(aov_LOH_Tx_CHROM)

xsq_LOH_CHROM_W_C <- chisq.test(LOHcounts_ID_CHROM_mean %>% filter(Tx_ID == "W") %>% pull(total_LOH),
                               LOHcounts_ID_CHROM_mean %>% filter(Tx_ID == "C") %>% pull(total_LOH))
xsq_LOH_CHROM_W_D <- chisq.test(LOHcounts_ID_CHROM_mean %>% filter(Tx_ID == "W") %>% pull(total_LOH),
                               LOHcounts_ID_CHROM_mean %>% filter(Tx_ID == "D") %>% pull(total_LOH))


# Permutation test for each chromosome and combo of strains
# Measures the difference in means, so normalized to n clones
all_CHROM_LOHrate_perm <- data.frame(NULL)
for(tx in 1:nrow(Tx_combo)){
  # tx = 1
  all_Chr_perm <- data.frame(Combo = with(Tx_combo[tx, ], paste(Tx_1, Tx_2, sep = "-")), 
                             CHROM = chrom_bound_BY$CHROM, 
                             obsvStat = 0, critVal = 0, pVal = 0, rejectNull = 0)
  
  for(c in seq_along(chrom_bound_BY$CHROM)) {
    # c = 14
    chr = chrom_bound_BY$CHROM[c]
    Chr_LOHcounts_ID <- LOHcounts_ID_CHROM %>% filter(CHROM == chr)
    Chr_perm <- perm_test(Chr_LOHcounts_ID, 
                          cat_var = "Tx_name", 
                          cat_names = with(Tx_combo[tx, ], c(Tx_1, Tx_2)), 
                          response_var = "n_LOH", 
                          alpha = 0.05,
                          n_perms = 10000, alt_hyp = "two-tailed",
                          rtrn = "all", include_matrix = F)
    all_Chr_perm[c, 3:6] <- unlist(Chr_perm[1:4])
  }
  all_CHROM_LOHrate_perm <- rbind(all_CHROM_LOHrate_perm, all_Chr_perm)
}

# all_CHROM_LOHrate_perm %>% arrange(desc(obsvStat)) %>% head()
# all_CHROM_LOHrate_perm %>% arrange(pVal) %>% head()

all_CHROM_LOHrate_perm_slim <- all_CHROM_LOHrate_perm %>% filter(Combo != "Cas9-Drive")
all_CHROM_LOHrate_perm_slim <- BHcorrection(all_CHROM_LOHrate_perm_slim)

CHROM_perm_sig <- all_CHROM_LOHrate_perm_slim %>% filter(pVal < FDR)
all_CHROM_LOHrate_perm_slim$sig_lab_p <- ifelse(all_CHROM_LOHrate_perm_slim$rejectNull == 1, "*", "")
all_CHROM_LOHrate_perm_slim$sig_lab_BH <- ifelse(all_CHROM_LOHrate_perm_slim$rejectNull_BH == 1, "**", "")
all_CHROM_LOHrate_perm_slim$sig_lab_both <- paste0(all_CHROM_LOHrate_perm_slim$sig_lab_p, 
                                                ifelse(all_CHROM_LOHrate_perm_slim$rejectNull_BH == 1, ", ", ""),
                                                all_CHROM_LOHrate_perm_slim$sig_lab_BH)

combo_cols <- colsplit(all_CHROM_LOHrate_perm_slim$Combo, "-", c("Tx1", "Tx2"))
all_CHROM_LOHrate_perm_slim <- cbind(all_CHROM_LOHrate_perm_slim, combo_cols)

# all_CHROMrate_perm_order %>% ggplot(aes(x = pVal)) + geom_histogram(binwidth = 0.05)

all_CHROM_LOHrate_perm_slim$Tx1 <- factor(all_CHROM_LOHrate_perm_slim$Tx1, 
                                          levels = c("WT", "Cas9", "Drive"))
all_CHROM_LOHrate_perm_slim$Tx2 <- factor(all_CHROM_LOHrate_perm_slim$Tx2, 
                                          levels = c("WT", "Cas9", "Drive"))
all_CHROM_LOHrate_perm_slim$Combo <- factor(all_CHROM_LOHrate_perm_slim$Combo)
sig_CHROMrate <- all_CHROM_LOHrate_perm_slim %>% filter(rejectNull_BH == 1) 

LOH_xCHROM_aov <- aov(n_LOH ~ Tx_name + CHROM, data = LOHcounts_ID_CHROM)
summary(LOH_xCHROM_aov)
LOH_xCHROM_Tuk <- TukeyHSD(LOH_xCHROM_aov, "Tx_name")
LOH_xCHROM_Tuk

###############################################################################
# Main Fig: Chromosome rates by Tx !!!!########################################
LOHcounts_ID_CHROM_mean$n_CHROM <- as.numeric(LOHcounts_ID_CHROM_mean$CHROM)
LOHcounts_ID_CHROM_mean$Tx_n <- as.numeric(LOHcounts_ID_CHROM_mean$Tx_ID)

point_diff <- 0.16

LOHrate_CHROM_line <- LOHcounts_ID_CHROM_mean %>% 
  ggplot() + 
  annotate(geom = "rect", xmin = LOHcounts_ID_CHROM_mean$n_CHROM[c(F, T)] - 0.5, 
           xmax = LOHcounts_ID_CHROM_mean$n_CHROM[c(F, T)] + 0.5, ymin = 0, ymax = Inf, 
           fill = "grey90", alpha = 0.1) +
  # geom_errorbar(aes(x = Tx_name, ymin = mean_rate, ymax = mean_rate, color = Tx_name), 
  #               width = 0.6, position = position_dodge(width = 0.6), size = 0.75) +
  geom_line(aes(x = n_CHROM + (Tx_n - 2) * point_diff, 
                y = bp_rate, color = Tx_ID, group = Tx_ID), 
             alpha = 0.7, size = 1) +
  geom_errorbar(aes(x = n_CHROM + (Tx_n - 2) * point_diff, 
                    ymin = bp_CI_lo, 
                    ymax = bp_CI_up, color = Tx_ID),
                width = 0, alpha = 0.7, size = 0.5) +
  geom_point(aes(x = n_CHROM + (Tx_n - 2) * point_diff, y = bp_rate, color = Tx_ID), 
             size = 4) +
  # geom_bracket(data = sig_CHROMrate,
  #              y.position = y_loc*1.05 + (as.numeric(sig_CHROMrate$Combo)-1) * y_loc * 0.03, step.increase = 0,
  #              vjust = 0.4, label.size = 5,
  #              aes(xmin = Tx1, xmax = Tx2, label = sig_lab_BH)) +
  # geom_text(aes(x = Tx_name, y = -0.00005, label = round(mean_LOH, 2), color = Tx_name),
  #           position = position_dodge(width = 0.6), show.legend = FALSE) +
  scale_color_manual(values = txPal, name = "Strain") +
  scale_y_continuous(labels = c(0, format(seq(5e-10, 2.5e-9, 5e-10), scientific = T)),
                     breaks = c(0, seq(5e-10, 2.5e-9, 5e-10))) +
  scale_x_continuous(labels = roman_chr, expand = c(0, 0),
                     breaks = 1:16, limits = c(0.5, 16.5)) +
  # ylim(c(0, NA)) +
  ylab("LOH event rate (/bp/gen)") + 
  xlab("Chromosome") +
  theme(axis.text.x = element_text(vjust = 3),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        # panel.background = element_rect(color = "grey80"),
        text = element_text(size = 24),
        legend.position = c(0.9, 0.9),
        legend.background = element_rect(fill = "white", color = "white"))

LOHrate_CHROM_line

ggsave(file.path(outIntDir, "LOHrate_CHROM_line_2022_04.png"), 
       plot = LOHrate_CHROM_line,
       device = "png",
       width = 16, height = 9, 
       units = "in",
       dpi = 600)

###############################################################################
# Empirical cumulative distribution of LOH rates ##############################
# along the length of each chromosome with permuted Anderson-Darling test
POSi_test_in <- POSi_data_in # %>%

# Hegazy test for all chromosomes and treatments against uniform !!!!#####
all_CHROM_uni_KS <- data.frame(NULL)
for(tx in 1:length(Tx_ID_levels)){
  # tx = 1
  tx_id <- Tx_ID_levels[tx]
  Combo_Chr_KS <- data.frame(Tx = tx_id,
                             CHROM = chrom_bound_BY$CHROM, 
                             KS_stat = 0, p_value = 0)
  KS_LOHcounts <- POSi_test_in %>% 
    filter(Tx_ID == tx_id) %>%
    select(Tx_ID, ID, CHROM, est_mid_POS)
  KS_LOHcounts$Tx_ID <- droplevels(KS_LOHcounts$Tx_ID)
  for(ch in seq_along(chrom_bound_BY$CHROM)) {
    # ch = 2
    chr = chrom_bound_BY$CHROM[ch]
    Chr_LOHcounts_ID <- KS_LOHcounts %>% filter(CHROM == chr)
    Chr_KS <- hegazy.unif.test(Chr_LOHcounts_ID$est_mid_POS/chrom_lengths_BY[ch])
    Combo_Chr_KS[ch, "KS_stat"] <- Chr_KS$statistic
    Combo_Chr_KS[ch, "p_value"] <- Chr_KS$p.value
  }
  all_CHROM_uni_KS <- rbind(all_CHROM_uni_KS, Combo_Chr_KS)
}

all_CHROM_uni_KS <- BHcorrection(all_CHROM_uni_KS, p_col = "p_value")

x <- c(0.33, 0.66, 1)

# Wasserstein test for all chromosomes and treatment combos !!!!#####
all_CHROM_LOHprop_AD <- data.frame(NULL)
for(tx in 1:nrow(Tx_ID_combo)){
  # tx = 1
  Tx_1 <- Tx_ID_combo$Tx_1[tx]
  Tx_2 <- Tx_ID_combo$Tx_2[tx]
  Combo_Chr_AD <- data.frame(Combo = with(Tx_ID_combo[tx, ], paste(Tx_1, Tx_2, sep = "-")),
                             Tx_1 = Tx_1, Tx_2 = Tx_2,
                             CHROM = chrom_bound_BY$CHROM, 
                             ADstat = 0, pVal = 0)
  LOH_1 <- POSi_test_in %>% 
    filter(Tx_ID == Tx_1) %>% 
    filter(!isTerm) %>%
    select(Tx_ID, ID, CHROM, est_mid_POS)
  LOH_2 <- POSi_test_in %>% 
    filter(Tx_ID == Tx_2) %>% 
    filter(!isTerm) %>%
    select(Tx_ID, ID, CHROM, est_mid_POS)
  # AD_LOHcounts$Tx_ID <- droplevels(AD_LOHcounts$Tx_ID)
  for(ch in seq_along(chrom_bound_BY$CHROM)) {
    # ch = 8
    chr = chrom_bound_BY$CHROM[ch]
    Chr_LOH_1 <- LOH_1 %>% filter(CHROM == chr)
    Chr_LOH_2 <- LOH_2 %>% filter(CHROM == chr)
    Chr_AD <- wass_test(Chr_LOH_1$est_mid_POS, Chr_LOH_2$est_mid_POS, nboots = 5000)
    Combo_Chr_AD[ch, c(5, 6)] <- Chr_AD[1:2]
  }
  all_CHROM_LOHprop_AD <- rbind(all_CHROM_LOHprop_AD, Combo_Chr_AD)
}

all_CHROM_LOHprop_AD_slim <- all_CHROM_LOHprop_AD %>% filter(Combo != "C-D")
all_CHROM_LOHprop_AD_slim <- BHcorrection(all_CHROM_LOHprop_AD_slim)
all_CHROM_LOHprop_AD_slim$p_adjust <- p.adjust(all_CHROM_LOHprop_AD_slim$pVal, method = "BH")

ad_label <- all_CHROM_LOHprop_AD_slim %>% filter(rejectNull_BH == 1)
ad_label$sig <- "*"
ad_label$Tx_2 <- factor(ad_label$Tx_2, levels = Tx_name_levels)
ad_label$rom_CHROM <- factor(chrom_IDs$rom_CHROM[as.numeric(ad_label$CHROM)],
                             levels = levels(chrom_IDs$rom_CHROM))

nLOH_CHROM <- POSi_test_in %>% count(Tx_name, CHROM)

marker_set <- data.frame(POSi = clean_markers)
marker_set_POS <- ConvertPosIndicies(marker_set, pos_col = "POSi", 
                                     index_out = "POS", add_chroms = T)
marker_set <- cbind(marker_set, marker_set_POS)

POSi_test_in$rom_CHROM <- factor(chrom_IDs$rom_CHROM[as.numeric(POSi_test_in$CHROM)],
                                 levels = levels(chrom_IDs$rom_CHROM))

# ad_label gives chromosomes and treatment pairs that are significant
# get position information to make brackets for indicating significance on plot

ecd_sig_region <- POSi_test_in %>%
  filter(!isTerm, rom_CHROM %in% ad_label$rom_CHROM, Tx_name %in% c("WT", "Drive")) %>% 
  arrange(Tx_name, est_mid_POS) %>% group_by(Tx_name) %>% 
  summarize(rom_CHROM = rom_CHROM, POS = est_mid_POS, f_bp = row_number()/n()) 

# Chr V is significant b/t WT and Drive. Cannot program bracket values
# for more than one significant combination, so this is coded for only one.
# find position that is maximally different between strains
ecd_sig_region <- ecd_sig_region %>% ungroup() %>% arrange(POS) %>% 
  mutate(lag_diff = f_bp - lag(f_bp))

ecd_sig_hi_POS <- ecd_sig_region %>% arrange(desc(lag_diff)) %>% head(1)
max_POS <- ecd_sig_hi_POS$POS + 5000

ecd_sig_region_lo <- ecd_sig_region %>% filter(Tx_name != ecd_sig_hi_POS$Tx_name)
# Find corresponding POS for lower value treatment
ecd_sig_lo_POS <- ecd_sig_region_lo[which.min(abs(max_POS - ecd_sig_region_lo$POS)), ]
ecd_bracket <- data.frame(rom_CHROM = ecd_sig_hi_POS$rom_CHROM, 
                          POS = b_POS, y_lo = ecd_sig_lo_POS$f_bp, y_hi = ecd_sig_hi_POS$f_bp)
ecd_bracket <- rbind(ecd_sig_hi_POS, ecd_sig_lo_POS)
ecd_bracket$POS <- mean(ecd_bracket$POS)

LOH_CHROM_ecd <- POSi_test_in %>%
  ggplot() +
  geom_segment(data = chrom_lengths_BY_df, 
               aes(x = 0, y = 0, xend = chrom_length/1000, yend = 1), size = 0.1, alpha = 0.9) +
  # geom_polygon(data = ecd_sig_region, aes(x = POS/1000, y = n, fill = Tx_name)) +
  stat_ecdf(aes(x = est_mid_POS/1000, color = Tx_ID)) + 
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
  facet_wrap(~rom_CHROM, ncol = 4, scales = "free_x") +
  theme(
    panel.grid.minor.y = element_blank(), 
    panel.grid.minor.x = element_blank(),
    panel.background = element_rect(color="grey80"),
    text = element_text(size = 18), legend.position = "bottom")

LOH_CHROM_ecd

ggsave(file.path(outIntDir, "iLOH_BPrate_xCHROM_ecd_2022_02.png"), 
       plot = LOH_CHROM_ecd,
       device = "png",
       width = 16, height = 9, 
       units = "in",
       dpi = 600)

chrom_in <- c("02", "03", "05", "09")

ecd_in <- POSi_data_in %>% filter(CHROM %in% chrom_in)
ecd_in$CHROM <- factor(ecd_in$CHROM, levels = chrom_in)
marker_in <- marker_set %>% filter(CHROM %in% chrom_in)
marker_in$CHROM <- factor(marker_in$CHROM, levels = chrom_in)
chrom_lengths_in <- chrom_lengths_BY_df %>% filter(CHROM %in% chrom_in)
chrom_lengths_in$CHROM <- factor(chrom_lengths_in$CHROM, levels = chrom_in)

LOH_CHROM_ecd_sub <- ecd_in %>%
  ggplot() +
   # geom_label(aes(x = chrom_bound_BY$End[1], y = 0.85, label = ad_label), hjust = "left") +
  # facet_wrap(~CHROM, ncol = 2, scales = "free_x") +
  stat_ecdf(aes(x = est_mid_POS/1000, color = Tx_name)) + 
  geom_point(data = marker_in, aes(x = POS/1000, y = -0.05), shape = "|", size = 1) +
  geom_segment(data = chrom_lengths_in, 
               aes(x = 0, y = 0, xend = chrom_length/1000, yend = 1), 
               size = 0.1, alpha = 0.9) +
  # stat_ecdf(aes(x = POSi), color = "red", size = 0.3) + 
  scale_color_manual(values = txPal, name = "Strain") + 
  xlab("Chromosome Position (kb)") + ylab("Cumulative Fraction of Breakpoints") +
  scale_x_continuous(expand = c(0.01, 0.01)) +
  facet_wrap(~CHROM, ncol = 2, scales = "free_x") +
  theme(
        panel.grid.minor.y = element_blank(), 
        panel.grid.minor.x = element_blank(),
        panel.background = element_rect(color="grey80"),
        text = element_text(size = 16), legend.position = "bottom")

LOH_CHROM_ecd_sub

ggsave(file.path(outIntDir, "iLOH_BPrate_ecd_sub_2022_02.png"), 
       plot = LOH_CHROM_ecd_sub,
       device = "png",
       width = 11, height = 8.5, 
       units = "in",
       dpi = 600)


###############################################################################
# LOH event rate for 1/3 chromosomes ####

all_LOHrates_Chrm3 <- POSi_data_wHet %>% 
  LOHrateSW(by_GT = F, win = "CHROM/3", chrom_df = chrom_bound_BY)

low_cover_ID_CHROM3 <- all_LOHrates_Chrm3 %>% 
  filter(prop_valid < 0.4) %>% select(ID, CHROM3) %>% 
  mutate(ID_CHROM3 = paste0(ID, "_", CHROM3))

chrom3s <- floor(chrom_lengths_BY/3)

chrom3_bound_BY <- data.frame(CHROM = chrom_bound_BY$CHROM, 
                              Start = chrom_bound_BY$Start,
                              mid1 = chrom_bound_BY$Start + chrom3s,
                              mid2 = chrom_bound_BY$Start + chrom3s*2,
                              End = chrom_bound_BY$End)

chrom3_len <- chrom3_bound_BY %>% 
  mutate(d_1 = mid1 - Start, d_2 = mid2 - mid1, d_3 = End - mid2) %>%
  select(!Start:End) %>%
  pivot_longer(cols = d_1:d_3, names_to = "chrom3", values_to = "len") %>%
  mutate(CHROM3 = paste0(CHROM, "_", substr(chrom3, 3, 3)))

POSi_test_in$CHROM3 <- "00"
for(i in 1:16) {
  # i = 1
  chr <- chrom3_bound_BY[i, "CHROM"]
  chr3 <- paste(chr, 1:3, sep = "_")
  i_pos1 <- POSi_test_in$est_mid >= chrom3_bound_BY[i, "Start"] & 
    POSi_test_in$est_mid < chrom3_bound_BY[i, "mid1"]
  i_pos2 <- POSi_test_in$est_mid >= chrom3_bound_BY[i, "mid1"] & 
    POSi_test_in$est_mid < chrom3_bound_BY[i, "mid2"]
  i_pos3 <- POSi_test_in$est_mid >= chrom3_bound_BY[i, "mid2"] & 
    POSi_test_in$est_mid <= chrom3_bound_BY[i, "End"]
  POSi_test_in$CHROM3[i_pos1] <- chr3[1]
  POSi_test_in$CHROM3[i_pos2] <- chr3[2]
  POSi_test_in$CHROM3[i_pos3] <- chr3[3]
}

POSi_test_in$CHROM3 <- factor(POSi_test_in$CHROM3)
POSi_test_in$ID_CHROM3 <- paste0(POSi_test_in$ID, "_", POSi_test_in$CHROM3)

POSi_test_in$lc_CHROM3 <- POSi_test_in$ID_CHROM3 %in% low_cover_ID_CHROM3$ID_CHROM3

LOHcounts_ID_CHROM3 <- POSi_test_in %>% 
  filter(!lc_CHROM3) %>%
  group_by(ID, CHROM3, .drop = F) %>% count(name = "n_LOH") %>% as.data.frame()
LOHcounts_ID_CHROM3 <- CategoriesFromID(LOHcounts_ID_CHROM3)

LOHcounts_ID_CHROM3_mean <- LOHcounts_ID_CHROM3 %>% 
  group_by(Tx_ID, CHROM3) %>% 
  summarize(total_LOH = sum(n_LOH), mean_LOH = mean(n_LOH), se_LOH = se(n_LOH))

LOHcounts_ID_CHROM3_mean$CHROM <- factor(substr(LOHcounts_ID_CHROM3_mean$CHROM3, 1, 2))
LOHcounts_ID_CHROM3_mean$CHROM3_pos <- factor(substr(LOHcounts_ID_CHROM3_mean$CHROM3, 4, 4))
LOHcounts_ID_CHROM3_mean$CHROM3_len <- rep(chrom3_len$len, 3)

LOH_CHROM3_CI_list <- LOHcounts_ID_CHROM3 %>% 
  group_by(Tx_ID, CHROM3) %>%
  group_modify(~ data.frame(boot.ci(boot(.x$n_LOH,
                                         statistic = sum.fun, R = 1000), 
                                    type = "perc")$percent[4:5] %>% t())) %>%
  rename(CI_95lo = X1, CI_95up = X2) %>% 
  select(CI_95lo, CI_95up)

LOHcounts_ID_CHROM3_mean <- merge(LOHcounts_ID_CHROM3_mean, 
                                 LOH_CHROM3_CI_list,
                                 by = c("Tx_ID", "CHROM3"))

LOHcounts_ID_CHROM3_mean <- LOHcounts_ID_CHROM3_mean %>% 
  mutate(bp_rate = total_LOH/CHROM3_len,
         bp_CI_lo = CI_95lo/CHROM3_len,
         bp_CI_up = CI_95up/CHROM3_len)

LOHcounts_ID_CHROM3_mean$bp_rate <- operate_by_factor_match(
  n_clones_LOH_xTx,
  LOHcounts_ID_CHROM3_mean[, c("Tx_ID", "bp_rate")],
  .fun = function(x, y) y/x)
LOHcounts_ID_CHROM3_mean$bp_CI_lo <- operate_by_factor_match(
  n_clones_LOH_xTx,
  LOHcounts_ID_CHROM3_mean[, c("Tx_ID", "bp_CI_lo")],
  .fun = function(x, y) y/x)
LOHcounts_ID_CHROM3_mean$bp_CI_up <- operate_by_factor_match(
  n_clones_LOH_xTx,
  LOHcounts_ID_CHROM3_mean[, c("Tx_ID", "bp_CI_up")],
  .fun = function(x, y) y/x)

# LOHcounts_ID_CHROM3_mean <- LOHcounts_ID_CHROM3_mean %>% mutate(bp_rate = mean_LOH/CHROM3_len, 
#                                                                 CI_95lo = (mean_LOH - se_LOH)/CHROM3_len, 
#                                                                 CI_95up = (mean_LOH + se_LOH)/CHROM3_len)
#   
ch3_names <- sort(unique(LOHcounts_ID_CHROM3$CHROM3))
all_CHROM3_LOHrate_perm <- data.frame(NULL)
for(tx in 1:nrow(Tx_combo)){
  # tx = 1
  all_chr3_perm <- data.frame(Combo = with(Tx_combo[tx, ], paste(Tx_1, Tx_2, sep = "-")), 
                              CHROM3 = ch3_names, 
                             obsvStat = 0, critVal = 0, pVal = 0, rejectNull = 0)
  for(c_i in seq_along(ch3_names)) {
    # chr = 10
    chr3 = ch3_names[c_i]
    chr3_LOHcounts_ID <- LOHcounts_ID_CHROM3 %>% filter(CHROM3 == chr3)
    chr3_perm <- perm_test(chr3_LOHcounts_ID, 
                          cat_var = "Tx_name", 
                          cat_names = with(Tx_combo[tx, ], c(Tx_1, Tx_2)), 
                          response_var = "n_LOH", 
                          alpha = 0.05,
                          n_perms = 10000, alt_hyp = "two-tailed",
                          rtrn = "all", include_matrix = F)
    all_chr3_perm[c_i, 3:6] <- unlist(chr3_perm[1:4])
    print(paste(all_chr3_perm[1, "Combo"], chr3, sep = " "))
  }
  all_CHROM3_LOHrate_perm <- rbind(all_CHROM3_LOHrate_perm, all_chr3_perm)
}

i_false_sig <- all_CHROM3_LOHrate_perm$pVal >= 0.05 & all_CHROM3_LOHrate_perm$rejectNull == 1
all_CHROM3_LOHrate_perm$rejectNull[i_false_sig] <- 0

all_CHROM3_LOHrate_perm_slim <- all_CHROM3_LOHrate_perm %>% filter(Combo != "Cas9-Drive")
all_CHROM3_LOHrate_perm <- BHcorrection(all_CHROM3_LOHrate_perm)
all_CHROM3_LOHrate_perm_slim <- BHcorrection(all_CHROM3_LOHrate_perm_slim)
all_CHROM3_LOHrate_perm_slim$p_adjust <- p.adjust(all_CHROM3_LOHrate_perm_slim$pVal)

all_CHROM3_LOHrate_perm$sig_lab_p <- ifelse(all_CHROM3_LOHrate_perm$rejectNull == 1, "*", "")
all_CHROM3_LOHrate_perm$sig_lab_BH <- ifelse(all_CHROM3_LOHrate_perm$rejectNull_BH == 1, "**", "")
all_CHROM3_LOHrate_perm$sig_lab_both <- paste0(all_CHROM3_LOHrate_perm$sig_lab_p, 
                                              ifelse(all_CHROM3_LOHrate_perm$rejectNull_BH == 1, ", ", ""),
                                              all_CHROM3_LOHrate_perm$sig_lab_BH)

combo_cols <- colsplit(all_CHROM3_LOHrate_perm$Combo, "-", c("Tx1", "Tx2"))
all_CHROM3_LOHrate_perm <- cbind(all_CHROM3_LOHrate_perm, combo_cols)

all_CHROM3_LOHrate_perm$Tx1 <- factor(all_CHROM3_LOHrate_perm$Tx1, levels = c("WT", "Cas9", "Drive"))
all_CHROM3_LOHrate_perm$Tx2 <- factor(all_CHROM3_LOHrate_perm$Tx2, levels = c("WT", "Cas9", "Drive"))
all_CHROM3_LOHrate_perm$Combo <- factor(all_CHROM3_LOHrate_perm$Combo)
all_CHROM3_LOHrate_perm$CHROM <- factor(substr(all_CHROM3_LOHrate_perm$CHROM, 1, 2))
all_CHROM3_LOHrate_perm$CHROM3_pos <- factor(substr(all_CHROM3_LOHrate_perm$CHROM, 4, 4))
sig_CHROM3rate <- all_CHROM3_LOHrate_perm %>% filter(rejectNull == 1)

y_loc <- max(LOHcounts_ID_CHROM3_mean$mean_rate + LOHcounts_ID_CHROM3_mean$CI_rate)
as.n <- function(x) {
  y <- as.numeric(x)
  return(y)
}

LOHrate_CHROM3_line <- LOHcounts_ID_CHROM3_mean %>% 
  ggplot() + 
  facet_wrap(~CHROM, ncol = 8, strip.position = "bottom", scales = "free_x") +
  geom_errorbar(aes(x = CHROM3_pos, y = bp_rate, ymin = bp_rate, ymax = bp_rate, color = Tx_ID), 
                width = 0.6, position = position_dodge(width = 0.85), size = 0.75) +
  geom_errorbar(aes(x = CHROM3_pos, ymin = bp_CI_lo,
                    ymax = bp_CI_up, color = Tx_ID),
                width = 0, position = position_dodge(width = 0.85), size = 0.75) +
  # geom_bracket(data = sig_CHROM3rate,
  #              y.position = y_loc*1.05 + (as.numeric(sig_CHROM3rate$Combo)-1)*y_loc*0.07,
  #              step.increase = 0, # step.group.by = sig_CHROM3rate$CHROM,
  #              vjust = 0, label.size = 5,
  #              aes(xmin = as.n(CHROM3_pos) + 0.333*0.85*(as.n(Tx1)-2),
  #                  xmax = as.n(CHROM3_pos) + 0.333*0.85*(as.n(Tx2)-2),
  #                  label = sig_lab_both)) +
  geom_text(aes(x = CHROM3_pos, y = -5E-8, label = total_LOH, color = Tx_ID),
            position = position_dodge(width = 0.85), show.legend = FALSE, size = 3.5) +
  scale_color_manual(values = txPal, name = "Strain") +
  # ylim(c(0, NA)) +
  ylab("LOH event rate (/bp)") + 
  xlab("Chromosome") +
  theme(axis.text.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.background = element_rect(color = "grey80"),
        text = element_text(size = 20))

LOHrate_CHROM3_line

ggsave(file.path(outIntDir, "LOHrate_CHROM3_line_2021_12.png"), 
       plot = LOHrate_CHROM3_line,
       device = "png",
       width = 16, height = 9, 
       units = "in",
       dpi = 600)

# LOHcounts_ID_CHROM3_mean %>% 
#   ggplot(aes(x = chr3, y  = mean_rate, color = Tx_name)) + 
#   geom_errorbar(aes(ymin = ifelse(mean_rate - se_rate < 0, 0, mean_rate - se_rate), 
#                     ymax = mean_rate + se_rate, color = Tx_name)) +
#   geom_point() +
#   # geom_line(aes(group = Tx_name)) +
#   scale_color_manual(values = txPal) + 
#   scale_fill_manual(values = txPal) + 
#   facet_wrap(~CHROM, ncol = 4, scales = "free_x")


LOHcounts_ID_CHROM3_mean %>% ggplot() + geom_histogram(aes(x = log(mean_rate)), bins = 5)

LOHcounts_ID_CHROM3_mean$ord  <- cut(LOHcounts_ID_CHROM3_mean$mean_rate, 
             breaks = nrow(LOHcounts_ID_CHROM3_mean), 
             include.lowest=TRUE,
             labels = str_pad(1:nrow(LOHcounts_ID_CHROM3_mean), width = 3, pad = "0"), ordered=TRUE)

LOHcounts_ID_CHROM3_mean %>% group_by(Tx_name) %>% summarize(v = var(mean_rate))
LOHcounts_ID_CHROM3_mean %>% group_by(chr3) %>% summarize(v = var(mean_rate)) %>% plot()

chr3_anova <- aov(mean_rate ~ chr3, data = LOHcounts_ID_CHROM3_mean)
summary(chr3_anova)

chr3_Tx_anova <- aov(mean_rate ~ Tx_name, data = LOHcounts_ID_CHROM3_mean)
summary(chr3_Tx_anova)

## fit ordered logit model and store results 'm'
m <- kruskal.test(mean_rate ~ Tx_name, data = LOHcounts_ID_CHROM3_mean)
m2 <- kruskal.test(mean_rate ~ chr3, data = LOHcounts_ID_CHROM3_mean)

m3 <- pairwise.wilcox.test(LOHcounts_ID_CHROM3_mean$mean_rate, LOHcounts_ID_CHROM3_mean$Tx_name)
m4 <- pairwise.wilcox.test(LOHcounts_ID_CHROM3_mean$mean_rate, 
                           LOHcounts_ID_CHROM3_mean$chr3, p.adjust.method = "none")

## view a summary of the model
summary(m)

###############################################################################
# Rate of LOH conversion for 50kb windows across each chromosome ####

all_LOHrates_SW50k <- POSi_data_wHet %>% 
  LOHrateSW(by_GT = F, win = 50000, slide_dist = 50000, chrom_df = chrom_bound_BY)
mean_LOHrates50k <- all_LOHrates_SW50k %>% 
  mutate(POSi = start_POSi, POS = start_POS) %>%
  rateMeans_xTx(per_gen = T) 
mean_LOHrates50k <- mean_LOHrates50k %>% 
  filter(n_clones >= 20)
mean_LOHrates50k <- ChromosomeCoordinates(mean_LOHrates50k)

all_LOHrates_SW50k_I <- POSi_data_wHet %>% 
  filter(!isTerm) %>% 
  LOHrateSW(by_GT = F, win = 50000, chrom_df = chrom_bound_BY)
mean_LOHrates50k_I <- all_LOHrates_SW50k_I %>% 
  rateMeans_xTx(per_gen = T)
mean_LOHrates50k_I <- mean_LOHrates50k_I %>% 
  filter(n_clones >= 20)
mean_LOHrates50k_I <- ChromosomeCoordinates(mean_LOHrates50k_I)

all_LOHrates_SW50k_T <- POSi_data_wHet %>% 
  filter(GT == "0/1" | isTerm) %>% 
  LOHrateSW(by_GT = F, win = 50000, chrom_df = chrom_bound_BY)

mean_LOHrates50k_T <- all_LOHrates_SW50k_T %>%
  rateMeans_xTx(per_gen = T) %>% 
  filter(n_sites > 0 & !is.na(rate))
mean_LOHrates50k_T <- mean_LOHrates50k_T %>% 
  filter(n_clones >= 20)
mean_LOHrates50k_T <- ChromosomeCoordinates(mean_LOHrates50k_T)

smooth_sw_in <- mean_LOHrates50k # %>% filter(n_sites >= 15000)
smooth_sw_in <- mean_LOHrates50k_I
smooth_sw_in$per_gen_rate <- smooth_sw_in$rate
smooth_sw_in$per_gen_se <- smooth_sw_in$std_err

smooth_sw_in$end_POS <- ConvertPosIndicies(smooth_sw_in, pos_col = "end_POSi", index_out = "POS")

smooth_sw_in$plot_POS <- (smooth_sw_in$POS + smooth_sw_in$end_POS)/2

lims <- smooth_sw_in %>% 
  group_by(CHROM) %>% 
  summarise(xmin=min(POS), xmax=max(end_POS))

lims$breaks <- floor((lims$xmax - lims$xmin)/50000)

ylims <- c(0, NA)

SW_LOH_plot <- ggplot() + 
  geom_rect(data = chrom_arms,
            aes(xmin = 0, xmax = (abs(arm_1_cent) + arm_2_cent)/1000,
                ymin = -Inf, ymax = Inf),
            alpha = 0.2, fill = "white", color = "grey70") +
  geom_vline(data = chrom_arms,
             aes(xintercept = abs(arm_1_cent)/1000), size = 0.35, color = "blue4", alpha = 0.5) +
  # geom_point(data = centrom_df[iCh, ], aes(x = POS, y = 0), color = "blue4", shape = 17, size = 1.5) +
  geom_ribbon(data = smooth_sw_in, 
              aes(x = plot_POS/1000, 
                  y = per_gen_rate,
                  xmin = plot_POS/1000, xmax = plot_POS/1000, 
                  ymin = ifelse((per_gen_rate - per_gen_se) < 0, 0, (per_gen_rate - per_gen_se)), 
                  ymax = (per_gen_rate + per_gen_se),
                  fill = Tx_name),
              alpha = 0.2) +
  geom_line(data = smooth_sw_in,
            aes(x = plot_POS/1000, y = per_gen_rate, color = Tx_name),
            size = 0.35)  +
  geom_point(data = smooth_sw_in,
             aes(x = plot_POS/1000, y = per_gen_rate, color = Tx_name),
             size = 0.2)  +
  # geom_blank(data = subset(smooth_sw_in, CHROM == "04"), 
  #            aes(x = plot_POS, y = rate/n_gens, color = Tx_name)) +
  scale_color_manual(values = txPal, name = "Strain") +
  scale_fill_manual(values = txPal, name = "Strain") +
  scale_x_continuous(expand = expansion(0, 0)) +
  scale_y_continuous(expand = expansion(0, 5E-6), 
                     labels = function(x) format(x, scientific = TRUE),
                     limits = ylims
                     # breaks = c(0, 1E-4, 2E-4, 3E-4, 4E-4, 5E-4)
  ) +
  theme(panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(color = "grey40", size = 0.1),
        legend.position = "bottom") +
  xlab("Chromosome position (kbp)") + ylab("LOH rate (/bp/generation)") +
  facet_wrap(~CHROM, ncol = 4, strip.position = "top", scales = "free_x") 

SW_LOH_plot

ggsave(file.path(outIntDir, "bCHROM_LOHrate_plot_2021_11.png"), 
       plot = SW_LOH_plot,
       device = "png",
       width = 16, height = 9, 
       units = "in",
       dpi = 600)

###############################################################################
# LOH rate vs chromsome size ####

LOHcounts_ID_CHROM_mean$chrom_length <- rep(chrom_lengths_BY, 3)

# Expected chromosome rate given genomic rate and chromosome size
LOHcounts_CHROM <- LOHcounts_ID_CHROM %>% 
  group_by(CHROM) %>% 
  summarize(total = sum(n_LOH))
LOH_CHROM_CI_list <- LOHcounts_ID_CHROM %>% 
  group_by(CHROM) %>%
  group_modify(~ data.frame(boot.ci(boot(.x$n_LOH,
                                         statistic = sum.fun, R = 1000), 
                                    type = "perc")$percent[4:5] %>% t())) %>%
  rename(CI_95lo = X1, CI_95up = X2) %>% 
  select(CI_95lo, CI_95up)

LOHcounts_CHROM <- merge(LOHcounts_CHROM, 
                         LOH_CHROM_CI_list,
                         by = c("CHROM"))

LOHcounts_CHROM$chrom_length <- chrom_lengths_BY
LOHcounts_CHROM <- LOHcounts_CHROM %>% mutate(rate = total/chrom_length/sum(n_clones_LOH_xTx$n)/n_gens,
                                              rate_CI_95lo = CI_95lo/chrom_length/sum(n_clones_LOH_xTx$n)/n_gens,
                                              rate_CI_95up = CI_95up/chrom_length/sum(n_clones_LOH_xTx$n)/n_gens)


expected_CHROM_rate <- lapply(chrom_lengths_BY_df$chrom_length,
                              function(y) sum(mean_LOHrate$sum_LOH)*(y/g_length)) %>%
  do.call("rbind", .) %>%
  as.data.frame() %>% rename(n_LOH = V1)

LOHcounts_CHROM <- cbind(LOHcounts_CHROM, expected_LOH = expected_CHROM_rate$n_LOH)

max(LOHcounts_CHROM$rate)/min(LOHcounts_CHROM$rate)

expected_TxCHROM_rate <- LOHcounts_ID_CHROM_mean %>% 
  select(Tx_ID, CHROM, total_LOH, bp_rate, bp_CI_lo, bp_CI_up)

expected_TxCHROM_rate <- lapply(mean_LOHrate$sum_LOH, 
                                function(x) sapply(chrom_lengths_BY_df$chrom_length,
                                                   function(y) x * (y/g_length))) %>% 
  as.data.frame(col.names = Tx_ID_levels) %>% 
  mutate(CHROM = chrom_lengths_BY_df$CHROM) %>% 
  pivot_longer(cols = 1:3, names_to = "Tx_ID", values_to = "expected_LOH") %>%
  merge(expected_TxCHROM_rate, ., by = c("Tx_ID", "CHROM")) %>% 
  arrange(Tx_ID)

expected_TxCHROM_rate <- lapply(mean_LOHrate$sum_LOH, 
       function(x) sapply(chrom_lengths_BY_df$chrom_length,
                          function(y) x * (y/g_length) / y / n_gens /
                            mean_LOHrate$n_clones[mean_LOHrate$sum_LOH == x])) %>% 
  as.data.frame(col.names = Tx_ID_levels) %>% 
  mutate(CHROM = chrom_lengths_BY_df$CHROM) %>% 
  pivot_longer(cols = 1:3, names_to = "Tx_ID", values_to = "e_bp_rate") %>% 
  arrange(Tx_ID, CHROM) %>% select(e_bp_rate) %>% cbind(expected_TxCHROM_rate, .)

# expected_TxCHROM_rate$per_gen <- expected_TxCHROM_rate$per_clone/n_gens
expected_TxCHROM_rate$chrom_length <- rep(chrom_lengths_BY, 3)


lm_coeff_list <- LOHcounts_ID_CHROM_mean %>% 
  group_by(Tx_ID) %>% 
  group_map(~ lm(bp_rate ~ chrom_length, data = .x)$coefficients)

names(lm_coeff_list) <- Tx_ID_levels
lm_coeff_df <- do.call(rbind, lm_coeff_list) %>% as.data.frame()
colnames(lm_coeff_df) <- c("b", "m")
lm_coeff_df$Tx_ID <- rownames(lm_coeff_df)
lm_coeff_df <- lm_coeff_df %>% 
  mutate(.x = min(chrom_lengths_BY), 
         .xend = max(chrom_lengths_BY), 
         .y = min(chrom_lengths_BY) * m + b, 
         .yend = max(chrom_lengths_BY) * m + b)


LOHcounts_ID_CHROM_mean %>% 
  group_by(CHROM) %>%
  # summarize(.chrom_length = mean(chrom_length),
  #           .mean_LOH = mean(bp_rate)) %>%
  ggplot() +
  # geom_line(data = expected_TxCHROM_rate,
  #           aes(x = chrom_length, y = bp_rate, color = Tx_ID), size = 1) +
  geom_segment(data = lm_coeff_df, aes(x = .x, xend = .xend, y = .y, yend = .yend, color = Tx_ID), size = 1) +
  # stat_smooth(data = LOH_chrom_reg$model,
  #             aes_string(x = names(LOH_chrom_reg$model)[2], y = names(LOH_chrom_reg$model)[1]), 
  #             method = "lm", col = "grey20", size = 0.5) +
  geom_point(aes(x = chrom_length, y = bp_rate, color = Tx_ID), size = 3) +
  # geom_abline(aes(intercept = LOH_chrom_reg$coefficients[1], slope = LOH_chrom_reg$coefficients[2])) +

  # geom_text(check_overlap = T, hjust = 0, aes(x = 1.1E6, y = 0.15), 
  #           label = paste0("y = ", formatC(LOH_chrom_reg$coef[[2]], 2, format = "E"), 
  #                          "x + ", signif(LOH_chrom_reg$coef[[1]], 3),
  #                         "\nAdj R2 = ", signif(summary(LOH_chrom_reg)$adj.r.squared, 3),
  #                    # "\nIntercept =",signif(LOH_chrom_reg$coef[[1]], 3),
  #                    # "\nSlope =", formatC(LOH_chrom_reg$coef[[2]], 2, format = "E"),
  #                    "\np = ",formatC(summary(LOH_chrom_reg)$coef[2,4], 2, format = "E"))) +
  scale_color_manual(values = txPal) +
  ylim(0, NA) + xlim(0,NA) 

ggplotReg <- function (fit) {
  
  require(ggplot2)
  
  ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + 
    geom_point() +
    stat_smooth(method = "lm", col = "red") +
    labs(title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
                       "Intercept =",signif(fit$coef[[1]],5 ),
                       " Slope =",signif(fit$coef[[2]], 5),
                       " P =",signif(summary(fit)$coef[2,4], 5)))
}

ggplotReg(LOH_chrom_reg)


###############################################################################
# Rate of breakpoints per clone per window within chromosomes !!! #############

POSi_data_in <- POSi_data_in %>% mutate(bp = 1)
LOH_mid_SW <- POSi_data_in %>% 
  mutate(POSi = est_mid) %>%
  # filter(!isTerm) %>%
  # filter(!cross_h_term) %>% 
  SliderCalc(data_col = "bp", index_col = "POSi", 
             window_size = 50000, slide_interval = 50000,
             summary_stat = sum, factor_col = "Tx_ID", chrom_win = T)
LOH_mid_SW$Tx_ID <- factor(LOH_mid_SW$Tx_ID, levels = Tx_ID_levels)


LOH_ID_SW <- POSi_data_in %>% 
  mutate(POSi = est_mid) %>%
  # filter(!isTerm) %>%
  # filter(!cross_h_term) %>% 
  SliderCalc(data_col = "bp", index_col = "POSi", 
             window_size = 50000, slide_interval = 50000,
             summary_stat = sum, factor_col = "ID", chrom_win = T)

LOH_ID_SW$Tx <- factor(substr(LOH_ID_SW$ID, 1, 1))
LOH_ID_SW$Tx_ID <- Recode_Tx_ID(LOH_ID_SW$Tx)
LOH_ID_SW$Tx_name <- Recode_Tx(LOH_ID_SW$Tx)
LOH_ID_SW <- LOH_ID_SW %>% select(-Tx)

fract_chrom_fun <- function(df, reps = 10000, alpha = 0.05) {
  # df <- LOH_ID_SW
  t_lo <- floor(reps * (alpha/2))
  t_up <- ceiling(reps * (1 - (alpha/2)))
  CI_df <- data.frame(NULL)
  for(tx in Tx_name_levels) {
    # tx <- "WT"
    n_id <- n_clones_xTx %>% filter(Tx_name == tx) %>% pull(n)
    tx_CI_df <- data.frame(NULL)
    for(chr in chrom_IDs$CHROM) {
      # chr <- "02"
      sub_df <- df %>% filter(Tx_name == tx, CHROM == chr)
      s_POSi <- sub_df %>% distinct(start) %>%
        arrange(start) %>% pull(start)
      n_s <- length(s_POSi)
      fract_m <- matrix(nrow = n_s, ncol = reps)
      for(rep in 1:reps) {
        # rep <- 8
        sample_df <- sub_df %>% group_by(start) %>% sample_n(size = n_id, replace = T)
        n_loh <- sum(sample_df$bp)
        n_loh <- ifelse(n_loh == 0, 1, n_loh)
        rep_v <- sample_df %>% group_by(start) %>% summarize(f = sum(bp)/n_loh) %>% pull(f)
        fract_m[, rep] <- rep_v
        # colnames(rep_df)[2] <- paste0("r_", rep)
        # fract_df <- cbind(fract_df, rep_df)
        
      }
      fract_m_sort <- t(apply(fract_m, 1, sort))
      chrom_CI_df <- data.frame(Tx_name = tx, CHROM = chr, start = s_POSi, 
                          CI_lo = fract_m_sort[, t_lo], CI_up = fract_m_sort[, t_up])
      tx_CI_df <- rbind(tx_CI_df, chrom_CI_df)
    }
    CI_df <- rbind(CI_df, tx_CI_df)
  }
  CI_df$Tx_name <- factor(CI_df$Tx_name, levels = Tx_name_levels)
  CI_df <- CI_df %>% arrange(Tx_name, start)
  return(CI_df)
}

LOH_SW_CIs <- fract_chrom_fun(LOH_ID_SW)
iLOH_SW_CIs <- LOH_SW_CIs

# colnames(LOH_SW_CIs)[grep("start", colnames(LOH_SW_CIs))] <- "start"
LOH_CHROM_counts <- POSi_data_in %>% filter(GT != "0/1", !isTerm) %>% 
  group_by(Tx_name, CHROM) %>% summarize(tot_LOH = n())

for(tx in Tx_name_levels) {
  # tx = "WT"
  # x = sub_counts$CHROM[1]
  i_SW <- LOH_mid_SW$Tx_name == tx
  sub_df <- LOH_mid_SW %>% filter(Tx_name == tx)
  sub_counts <- LOH_CHROM_counts %>% filter(Tx_name == tx)
  fract_list <- sapply(sub_counts$CHROM, function(x) 
    sub_df$bp[sub_df$CHROM == x]/sub_counts$tot_LOH[sub_counts$CHROM == x])
  LOH_mid_SW$chrom_fraction[i_SW] <- unlist(fract_list)
}


# LOH_mid_SW <- LOH_mid_SW %>% group_by(Tx_name, CHROM) %>% mutate(chrom_fraction = bp/sum(bp))

# LOH_mid_SW <- merge(LOH_mid_SW, LOH_SW_CIs, by = c("Tx_name", "CHROM", "start")) %>% arrange(Tx_name, start)

# Per clone rate. Each Tx has a different # of clones
# LOH_mid_SW$clone_rate <- 0
# LOH_mid_SW$clone_rate[LOH_mid_SW$Tx_name == "WT"] <- 
#   LOH_mid_SW$bp[LOH_mid_SW$Tx_name == "WT"]/n_evo_xTx$n[n_evo_xTx$Tx_name == "WT"]
# LOH_mid_SW$clone_rate[LOH_mid_SW$Tx_name == "Cas9"] <- 
#   LOH_mid_SW$bp[LOH_mid_SW$Tx_name == "Cas9"]/n_evo_xTx$n[n_evo_xTx$Tx_name == "Cas9"]
# LOH_mid_SW$clone_rate[LOH_mid_SW$Tx_name == "Drive"] <- 
#   LOH_mid_SW$bp[LOH_mid_SW$Tx_name == "Drive"]/n_evo_xTx$n[n_evo_xTx$Tx_name == "Drive"]

LOH_mid_SW <- LOH_mid_SW %>% mutate(start_POSi = round(start), end_POSi = round(end)) %>% select(!start:end)
s_POS <- LOH_mid_SW %>% ConvertPosIndicies(pos_col = "start_POSi", index_out = "POS", add_chroms = T)
names(s_POS)[1] <- "start_POS"
e_POS <- LOH_mid_SW %>% ConvertPosIndicies(pos_col = "end_POSi", index_out = "POS")
LOH_mid_SW <- data.frame(LOH_mid_SW, s_POS[c(2,1)], end_POS = e_POS)
colnames(LOH_mid_SW)[1] <- c("Tx_name")
LOH_mid_SW$Tx_name <- factor(LOH_mid_SW$Tx_name, levels = c("WT", "Cas9", "Drive"))
LOH_mid_SW <- LOH_mid_SW %>% 
  mutate(mid_POS = start_POS + 50000/2, mid_POSi = (start_POSi + end_POSi)/2,
         mid_POS_kb = (start_POS + 50000/2)/1000, mid_POSi_kb = (start_POSi + end_POSi)/2000)

LOH_mid_SW$rom_CHROM <- factor(chrom_IDs$rom_CHROM[as.numeric(LOH_mid_SW$CHROM)], levels = roman_chr)
centrom_df$rom_CHROM <- factor(chrom_IDs$rom_CHROM[as.numeric(centrom_df$CHROM)], levels = roman_chr)

breakpoints_SW_plot <- LOH_mid_SW %>% 
  # filter(CHROM == "03") %>%
  # filter(Tx_ID == "D") %>%
  ggplot(aes(x = mid_POS_kb, y = chrom_fraction)) + 
  geom_vline(data = centrom_df, aes(xintercept = POS/1000), 
             linetype = 2,color = "blue4", size = 0.25) +
  # geom_ribbon(aes(ymin = CI_lo, ymax = CI_up, fill = Tx_name, color = Tx_name), 
  #             size = 0.05, alpha = 0.1) +
  geom_line(aes(color = Tx_ID)) +
  # geom_line(data = genome_scores_50kb, 
  #           aes(x = mid_POS/1000, y = -(score - min_score)/((max(score) - min_score)*3))) +
  # geom_line(aes(x = mid_POS_kb, y = log(rate, 10), color = Tx_name)) +
  # geom_point(aes(x = mid_POS_kb, y = rate, color = Tx_name)) +
  # xlim(0, max(chrom_lengths_BY)) +
  xlab("Window Position (kb)") + ylab("Breakpoints in window / Breakpoints in chromosome") +
  xlim(0,NA) +
  # scale_y_continuous(limits = c(-1/3, 0.4), 
  #                    breaks = c(-1/3, -2/9, -1/9, seq(0, 0.4, 0.1)), 
  #                    labels = c(210, 140, 70, seq(0, 0.4, 0.1))) +
  scale_color_manual(values = txPal, name = "Strain") + 
  scale_fill_manual(values = txPal, name = "Strain") + 
  facet_wrap(~rom_CHROM, ncol = 4, scales = "free_x") +
  theme(legend.text = element_text(size = 13),
        legend.title = element_text(size = 14),
        strip.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 16),
        # text = element_text(size = 18), 
        plot.margin = margin(20, 20, 20, 20),
        panel.grid.minor = element_blank(),
        legend.position = c(0.97, 0.94),
        legend.background = element_rect(fill = "white", color = "white"))
  
breakpoints_SW_plot

ggsave(file.path(outIntDir, "iLOH_bp_SW50x5_plot.png"), 
       plot = breakpoints_SW_plot,
       device = "png",
       width = 16, height = 9, 
       units = "in",
       dpi = 600)

POSi_data_in %>% 
  ggplot() + geom_histogram(aes(x = est_mid, fill = Tx_name), binwidth = 10000) +
  # geom_point(aes(x = mid_POS, y = value, color = Tx_name)) +
  # xlim(0, max(chrom_lengths_BY)) +
  xlab("Window Position (kbp)") + ylab("Number of Breakpoints") +
  xlim(0, NA) +
  scale_fill_manual(values = txPal, name = "Drive Type") + 
  facet_wrap(~CHROM, ncol = 4, scales = "free_x") +
  theme(text = element_text(size = 14), legend.position = "bottom")

###############################################################################
# Breakpoint /bp/clone rate vs position #####
LOH_BPrate_SW <- POSi_data_in %>%
  mutate(bp = 1) %>%
  # filter(!isTerm) %>%
  # filter(!cross_h_term) %>% 
  # mutate(POS = fract_dist) %>%
  mutate(POSi = est_mid) %>%
  SliderCalc(data_col = "bp", index_col = "POSi", 
             window_size = 50000, slide_interval = 5000,
             summary_stat = sum, factor_col = "Tx_ID", chrom_win = T)
LOH_BPrate_SW$Tx_ID <- factor(LOH_BPrate_SW$Tx_ID, levels = Tx_ID_levels)

LOH_BPrate_SW$start_POS <- LOH_BPrate_SW %>% ConvertPosIndicies(pos_col = "start", index_out = "POS")
LOH_BPrate_SW$end_POS <- LOH_BPrate_SW %>% ConvertPosIndicies(pos_col = "end", index_out = "POS")
LOH_BPrate_SW$mid_POS <- round((LOH_BPrate_SW$start_POS + LOH_BPrate_SW$end_POS)/2)
LOH_BPrate_SW$length <- LOH_BPrate_SW$end - LOH_BPrate_SW$start + 1
LOH_BPrate_SW$per_clone <- operate_by_factor_match(n_clones_LOH_xTx, 
                                                 LOH_BPrate_SW[, c("Tx_ID", "bp")],
                                                 .fun = function(x, y) y/x)
LOH_BPrate_SW <- LOH_BPrate_SW %>% 
  mutate(bp_rate = per_clone/length/n_gens)

LOH_BPrate_SW <- LOH_BPrate_SW %>% 
  mutate(mid_POSi = (start + end)/2,
         mid_POS_kb = mid_POS/1000, mid_POSi_kb = (start + end)/2000)

LOH_BPrate_SW$rom_CHROM <- factor(chrom_IDs$rom_CHROM[as.numeric(LOH_BPrate_SW$CHROM)], levels = roman_chr)

max_rate <- max(LOH_BPrate_SW$bp_rate) * 1.2
y_scale <- ceiling(log10(max_rate)) - 1
# y_max <- 10^y_scale
y_max <- 5*10^(y_scale-1)
y_incr <- 2.5*10^(y_scale - 1)

scale_bar <- data.frame(rom_CHROM = chrom_IDs$rom_CHROM, scale = 50, 
                        .x = 0, .xend = 50, .y = 5E-9, .yend = 5E-9)

BPrate_SW_plot <- LOH_BPrate_SW %>% 
  # filter(length > 10000) %>%
  # filter(bp < 10) %>%
  filter(bp_rate < 5E-9) %>%
  # filter(CHROM == "03") %>%
  # filter(Tx_ID == "D") %>%
  ggplot(aes(x = start_POS/1000, y = bp_rate)) + 
  geom_segment(data = scale_bar, 
               aes(x = .x, xend = .xend, y = .y, yend = .yend),
               size = 1, color = "grey30") +
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
  # xlab("Window Position (kb)") + xlim(0, NA) +
  scale_x_continuous(name = "Window Position (kb)", 
                     limits = c(0, NA)) +
  scale_y_continuous(labels = format_sci_10,
    # labels = c(0, format(seq(y_incr, y_max, y_incr), scientific = T)),
                     # breaks = c(0, seq(y_incr, y_max, y_incr)),
  name = "Breakpoint Rate /bp/gen") +
  scale_color_manual(values = txPal, name = "Strain") + 
  scale_fill_manual(values = txPal, name = "Strain") + 
  facet_wrap(~rom_CHROM, ncol = 4, scales = "free_x") +
  theme(legend.text = element_text(size = 16),
        legend.title = element_text(size = 20),
        strip.text = element_text(size = 20),
        axis.title.y = element_text(size = 24, vjust = 2.5),
        axis.title.x = element_text(size = 24, vjust = -1.7),
        axis.text = element_text(size = 18),
        # text = element_text(size = 18), 
        plot.margin = margin(20, 20, 20, 20),
        panel.grid.minor = element_blank(),
        legend.position = c(0.96, 0.96),
        legend.background = element_rect(fill = "white", color = "white"))

BPrate_SW_plot

ggsave(file.path(outIntDir, "BPrate_SW_SW50x5_plot_2022_05.png"), 
       plot = BPrate_SW_plot,
       device = "png",
       width = 16, height = 9, 
       units = "in",
       dpi = 600)

###############################################################################
# Breakpoint rate as a function of marker density in 50kb windows #############
LOH_BPrate_SW50_50 <- POSi_data_in %>%
  mutate(bp = 1) %>%
  # filter(!isTerm) %>%
  # filter(!cross_h_term) %>% 
  # mutate(POS = fract_dist) %>%
  mutate(POSi = est_mid) %>%
  SliderCalc(data_col = "bp", index_col = "POSi", 
             window_size = 50000, slide_interval = 50000,
             summary_stat = sum, factor_col = "Tx_ID", chrom_win = T)
LOH_BPrate_SW50_50$Tx_ID <- factor(LOH_BPrate_SW50_50$Tx_ID, levels = Tx_ID_levels)

LOH_BPrate_SW50_50$start_POS <- LOH_BPrate_SW50_50 %>% ConvertPosIndicies(pos_col = "start", index_out = "POS")
LOH_BPrate_SW50_50$end_POS <- LOH_BPrate_SW50_50 %>% ConvertPosIndicies(pos_col = "end", index_out = "POS")
LOH_BPrate_SW50_50$mid_POS <- round((LOH_BPrate_SW50_50$start_POS + LOH_BPrate_SW50_50$end_POS)/2)
LOH_BPrate_SW50_50$length <- LOH_BPrate_SW50_50$end - LOH_BPrate_SW50_50$start + 1
LOH_BPrate_SW50_50$per_clone <- operate_by_factor_match(n_clones_LOH_xTx, 
                                                   LOH_BPrate_SW50_50[, c("Tx_ID", "bp")],
                                                   .fun = function(x, y) y/x)
LOH_BPrate_SW50_50 <- LOH_BPrate_SW50_50 %>% 
  mutate(bp_rate = per_clone/length/n_gens)

LOH_BPrate_SW50_50 <- LOH_BPrate_SW50_50 %>% 
  mutate(mid_POSi = (start + end)/2,
         mid_POS_kb = mid_POS/1000, mid_POSi_kb = (start + end)/2000)

LOH_BPrate_SW50_50$rom_CHROM <- factor(chrom_IDs$rom_CHROM[as.numeric(LOH_BPrate_SW50_50$CHROM)], levels = roman_chr)

LOH_BPrate_SW50_50$Tx_ID <- Recode_Tx_ID(LOH_BPrate_SW50_50$Tx_ID, "Tx_ID")

marker_SW_Tx <- LOH_SNPs %>%
  mutate(m = 1) %>%
  SliderCalc(data_col = "m", index_col = "POSi", 
             window_size = 50000, slide_interval = 50000,
             summary_stat = sum, factor_col = "Tx_ID", chrom_win = T) 


marker_SW_Tx <- marker_SW_Tx %>% mutate(start_POSi = round(start), end_POSi = round(end), mid_POSi = (start + end)/2) %>% select(!start:end)
# s_POS <- marker_SW_Tx %>% ConvertPosIndicies(pos_col = "start_POSi", index_out = "POS", add_chroms = T)
# names(s_POS)[1] <- "start_POS"
# e_POS <- marker_SW_Tx %>% ConvertPosIndicies(pos_col = "end_POSi", index_out = "POS")
# marker_SW_Tx <- data.frame(marker_SW_Tx, s_POS[c(2,1)], end_POS = e_POS)
marker_SW_Tx$Tx_ID <- Recode_Tx_ID(marker_SW_Tx$Tx_ID, "Tx_ID")
# marker_SW_Tx <- marker_SW_Tx %>% 
#   mutate(mid_POS = start_POS + 50000/2, mid_POSi = (start_POSi + end_POSi)/2,
#          mid_POS_kb = (start_POS + 50000/2)/1000, mid_POSi_kb = (start_POSi + end_POSi)/2000)
# marker_SW_Tx$rom_CHROM <- factor(chrom_IDs$rom_CHROM[as.numeric(marker_SW_Tx$CHROM)], levels = roman_chr)

marker_SW_Tx <- marker_SW_Tx %>%
  mutate(pc_m = operate_by_factor_match(n_clones_LOH_xTx, 
                                        data.frame(Tx_ID, m),
                                        function(x, y) y/x))

marker_SW_Tx <- merge(marker_SW_Tx[, c("Tx_ID", "mid_POSi", "m", "pc_m")], 
                      LOH_BPrate_SW50_50, 
                      by = c("Tx_ID", "mid_POSi")) %>%
  arrange(Tx_ID, mid_POSi)

marker_SW_Tx %>% ggplot() + geom_point(aes(x = pc_m, y = per_clone))

lm_m <- lm(bp ~ log10(m), data = marker_SW_Tx %>% filter(m > 0))
summary(lm_m)

marker_SW_Tx %>% ggplot(aes(x = m)) + 
  geom_histogram(aes(fill = Tx_ID), 
                 position = position_dodge()) +
  scale_fill_manual(values = txPal)

marker_SW_Tx %>% ggplot(aes(x = log10(m), y = bp)) + 
  geom_point(aes(color = Tx_ID)) +
  geom_smooth(method = "lm") +
  scale_x_continuous(limits = c(2, 5)) +
  scale_color_manual(values = txPal)


###########################################################################
# Breakpoint distribution against fractional distance from centromere #####

LOH_cent_SW <- POSi_data_in %>% 
  mutate(bp = 1) %>%
  # filter(!isTerm) %>%
  # filter(!cross_h_term) %>% 
  # mutate(POS = dist_cent_mid_abs) %>%
  mutate(POSi = dist_cent_mid_abs) %>%
  SliderCalc(data_col = "bp", index_col = "POSi", 
             window_size = 50000, slide_interval = 50000,
             summary_stat = sum, factor_col = "type", chrom_win = F)
LOH_cent_SW$Tx_ID <- factor(LOH_cent_SW$Tx_ID, levels = Tx_ID_levels)

LOH_cent_SW$start_POS <- LOH_cent_SW %>% ConvertPosIndicies(pos_col = "start", index_out = "POS")
LOH_cent_SW$end_POS <- LOH_cent_SW %>% ConvertPosIndicies(pos_col = "end", index_out = "POS")
LOH_cent_SW$mid_POS <- round((LOH_cent_SW$start_POS + LOH_cent_SW$end_POS)/2)
LOH_cent_SW$length <- LOH_cent_SW$end - LOH_cent_SW$start
LOH_cent_SW$per_clone <- operate_by_factor_match(n_clones_LOH_xTx, 
                                                 LOH_cent_SW[, c("Tx_ID", "bp")],
                                                 .fun = function(x, y) y/x)

# iLOH_CIs <- data.frame(NULL)
# for(tx in Tx_ID_levels) {
#   # tx <- "W"
#   tx_n_LOH <- LOH_cent_SW %>% filter(Tx_ID == tx) %>% pull(n_iLOH)
#   
#   tx_boot <- boot(tx_n_LOH, statistic = mean.fun, R = 10000)
#   tx_boot_CI <- boot.ci(boot.out = tx_boot, conf = 0.95,
#                         type = "basic")$basic[4:5]
#   tx_CI <- data.frame(Tx_ID = tx, low95CI = tx_boot_CI[1], up95CI = tx_boot_CI[2])
#   iLOH_CIs <- rbind(iLOH_CIs, tx_CI)
# }
# LOH_cent_SW <- merge(LOH_cent_SW, iLOH_CIs, by = "Tx_ID", all = T)

# df to normalize number of sites at each position sampled
left_arm <- centrom_df$POS
right_arm <- chrom_lengths_BY_df$chrom_length - centrom_df$POS
chrom_arms_df <- data.frame(POS = sort(c(left_arm, right_arm), decreasing = T), i = 1)
chrom_arms_df$n_sampled <- cumsum(chrom_arms_df$i)
LOH_cent_SW$bp_adj <- LOH_cent_SW$bp
for(i in 1:(nrow(chrom_arms_df) - 1)) {
  # i = 1
  i_adj <- LOH_cent_SW$mid_POS <= chrom_arms_df$POS[i]
  # LOH_cent_SW$n_sampled[i] <- LOH_cent_SW$bp[i]/chrom_arms_df$n_sampled[i]
  LOH_cent_SW$bp_adj[i_adj] <- LOH_cent_SW$bp[i_adj]/chrom_arms_df$n_sampled[i]
}


# Normalize to the number of bp on chromosome
LOH_cent_SW %>% 
  # filter(per_clone < 0.2) %>%
  ggplot(aes(x = start, y = bp_adj, color = type, group = type)) +
  geom_line() +
  geom_point(size = 0.25) + 
  # facet_wrap(~CHROM, ncol = 4, scales = "free") +
  # ylim(0, 0.0006) +
  scale_color_manual(values = txPal)

POSi_data_in %>% 
  ggplot(aes(x = dist_cent_mid, fill = type)) + 
  geom_histogram(binwidth = 50000, position = "stack") + 
  facet_wrap(~CHROM, ncol = 4, scales = "free_x")

POSi_data_in %>% 
  # filter(CHROM != "03", CHROM != "09") %>%
  ggplot(aes(x = fract_dist, fill = CHROM)) + 
  geom_histogram(binwidth = 0.05, position = "stack") + 
  facet_wrap(~type, ncol = 1) +
  scale_fill_manual(values = lineagePal)

POSi_data_in %>%
  group_by(CHROM, type) %>%
  summarize(mean_cent = mean(dist_cent_mid_abs),
            mean_term = mean(dist_term)) %>%
  as.data.frame()

# POSi_data_in %>% filter(CHROM == "09", est_mid > 5820000)
# POSi_data_in %>% filter(CHROM == "09", est_end_POS > 400000) %>% mutate(isTerm, in_POS = 439888 - est_end_POS)
# POSi_data_wHet %>% filter(CHROM == "09", est_end_POS > 400000) %>% head()
# 
# POSi_data_wHet %>% filter(CHROM == "09", est_end > 5760000) %>% head()
# POSi_data_wHet %>% filter(CHROM == "09", ID == "F_A02")

###############################################################################
# Breakpoint rate vs gRNA similarity in Drive !!! #############################
DrivevsWT_breaks <- LOH_mid_SW %>% 
  select(Tx_name, CHROM, rate, mid_POS_kb, mid_POSi_kb) %>%
  filter(Tx_name != "Cas9") %>% 
  pivot_wider(names_from = Tx_name, values_from = rate) %>%
  mutate(rate_diff = (Drive - WT)/(Drive + WT + 0.00001))

genome_scores_50kb <- genome_scores_50kb %>%
  mutate(mid_POS_kb = mid_POS/1000, norm_score = (score - min_score)/(max(score) - min_score))

DrivevsWT_POS_plot <- DrivevsWT_breaks %>% 
  # filter(CHROM == "05") %>%
  ggplot() + 
  geom_line(aes(x = mid_POS_kb, y = rate_diff)) +
  geom_line(data = genome_scores_50kb,
            aes(x = mid_POS_kb, y = norm_score), color = "red3") +
  # geom_line(aes(x = mid_POS_kb, y = log(rate, 10), color = Tx_name)) +
  # geom_point(aes(x = mid_POS_kb, y = rate, color = Tx_name)) +
  # xlim(0, max(chrom_lengths_BY)) +
  xlab("Window Position (kbp)") + ylab("Number of Breakpoints / Clone") +
  xlim(0,NA) +
  # scale_y_continuous(limits = c(-1/3, 0.3), 
  #                    breaks = c(-1/3, -2/9, -1/9, seq(0, 0.3, 0.1)), 
  #                    labels = c(210, 140, 70, seq(0, 0.3, 0.1))) +
  # scale_color_manual(values = txPal, name = "Drive Type") + 
  facet_wrap(~CHROM, ncol = 4, scales = "free_x") +
  theme(text = element_text(size = 14), legend.position = "bottom")

DrivevsWT_POS_plot

DrivevsWT_breaks_scores <- cbind(DrivevsWT_breaks, norm_score = genome_scores_50kb$norm_score)

DrivevsWT_breaks_scores %>% 
  # filter(!(WT == 0 & Drive == 0)) %>%
  ggplot() + 
  geom_point(aes(x = rate_diff, y = norm_score)) +
  xlab("Breakpoint rate bias") + ylab("Normailzed gRNA similarity score") +
  ylim(0,1) +
  scale_x_continuous(breaks = c(-1, -0.5, 0, 0.5, 1), 
                     labels = c("WT only", "WT higher", "0", "Drive higher", "Drive only")) + 
  theme(panel.grid.minor.x = element_blank())


  
    