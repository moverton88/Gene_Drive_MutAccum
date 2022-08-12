
# Table of key overall rates
# n clones
# n LOH
# n true LOH
# true LOH rate /gen/gen (CI)
# n LOH to BY
# n true LOH to BY
# n true iLOH
# n true tLOH
# n PMs
# PM rate /bp/gen (CI)
# n indels
# indel rate /bp/gen (CI)

larger_CI <- function(df) {
  ci <- rep(0, nrow(df))
  if(length(df) > 2) {
    i_CI_lo <- df[, 1] - df[, 2] > df[, 3] - df[, 1]
    ci[i_CI_lo] <- (df[, 1] - df[, 2])[i_CI_lo]
    ci[!i_CI_lo] <- (df[, 3] - df[, 1])[!i_CI_lo]
  }
  return(ci)
}

get_E <- function(nums) {
  floor(log10(nums))
}

get_sci_digits <- function(nums, base = "data", n_digits = 2) {
  if(base == "data") {
    round(nums * 10 ^ -floor(log10(nums)), n_digits)
  } else {
    round(nums * 10 ^ -base, n_digits)
  }
}


tx_pLOH <- all_GT_bounds_merge %>% filter(GT != "0/1") %>% count(Tx_ID, name = "putative_LOH")
tx_pLOH_BY <- all_GT_bounds_merge %>% filter(GT == "0/0") %>% count(Tx_ID, name = "putative_BY")
tx_LOHrate <- tx_LOHrate %>% rename(n_clones_LOH = n_clones)
tx_SNMrate <- tx_SNMrate %>% rename(n_clones_SNM = n_clones)


tx_master_table <- cbind(tx_LOHrate[, 1:2],
                         tx_pLOH[2], tx_pLOH_BY[2],
                         tx_LOHrate[, 3:length(tx_LOHrate)],
                         tx_SNMrate[, 2:length(tx_SNMrate)],
                         tx_indelRate[, 3:length(tx_indelRate)])

tx_master_table <- tx_master_table %>% 
  mutate(BP_rate_sci = get_sci_digits(BP_rate),
         SNM_rate_sci = get_sci_digits(SNM_rate),
         indel_rate_sci = get_sci_digits(indel_rate),
         BP_E = get_E(BP_rate),
         SNM_E = get_E(SNM_rate),
         indel_E = get_E(indel_rate))

tx_master_table$BP_CI_sci <- round(larger_CI(tx_master_table %>% 
                                          select(contains("BP_"))) *
                               10^-tx_master_table$BP_E, 2)
tx_master_table$SNM_CI_sci <- round(larger_CI(tx_master_table %>% 
                                           select(contains("SNM_"))) *
                                10^-tx_master_table$SNM_E, 2)
tx_master_table$indel_CI_sci <- round(larger_CI(tx_master_table %>% 
                                             select(contains("indel_"))) *
                                  10^-tx_master_table$indel_E, 2)
tx_master_table$fLOH_CI_sci <- round(larger_CI(tx_master_table %>% 
                                                  select(contains("f_"))), 3)

tx_master_table$Strain <- Recode_Tx_ID(tx_master_table$Strain, "Tx_ID")
tx_master_table <- tx_master_table %>% arrange(Strain)

tx_master_in <- tx_master_table %>% 
  mutate(iLOH_label = paste0(total_iLOH, " (", round(100*f_iLOH, 1), "%)"),
         BY_label = paste0(total_BY, " (", round(100*(total_BY/total_LOH), 1), "%)")) %>%
  # rename(Strain = Tx_ID) %>%
  select(Strain, n_clones_LOH,
         total_LOH, BP_rate_sci,
         BY_label,
         # total_BY,  
         iLOH_label,
         # total_iLOH, 
         # f_iLOH, 
         # fLOH_CI_sci,
         n_clones_SNM,
         total_SNM, SNM_rate_sci, 
         total_indels, indel_rate_sci,
         BP_CI_sci:indel_CI_sci,
         contains("E", ignore.case = F))


master_table <- tx_master_in %>% 
  # select(!contains("CIlo") & !contains("CIup")) %>%
  select(!contains("putative") & !contains("tLOH")) %>%
  gt(rowname_col = "TX_ID")  %>% 
  # fmt_scientific(columns = matches("rate|CI"), decimals = 2) %>% 
  cols_merge_uncert(
    col_val = BP_rate_sci,
    col_uncert = BP_CI_sci) %>%
  # cols_merge_uncert(
  #   col_val = f_iLOH,
  #   col_uncert = fLOH_CI_sci) %>%
  cols_merge_uncert(
    col_val = SNM_rate_sci,
    col_uncert = SNM_CI_sci) %>%
  cols_merge_uncert(
    col_val = indel_rate_sci,
    col_uncert = indel_CI_sci) %>%
  fmt_number(columns = matches("_rate_"), decimals = 2) %>%
  # fmt_number(columns = matches("f_"), decimals = 3) %>%
  cols_label(n_clones_LOH = "End-point clones",
             # putative_LOH = "Total putative",
             # putative_BY = "Putative convert BY",
             total_LOH = "Total Events",
             BP_rate_sci = html("Event rate (&times10<sup>&#8722<span>10</sup>)"),
             BY_label = "RM/BY to BY/BY",
             # total_BY = "RM/BY to BY/BY",
             iLOH_label = "Total interstitial",
             # total_iLOH = "Total interstitial",
             # total_tLOH = "Total terminal",
             # f_iLOH = "Fraction interstitial",
             n_clones_SNM = "End-point clones",
             total_SNM = "Total",
             SNM_rate_sci = html("Rate (&times10<sup>&#8722<span>10</sup>)"),
             total_indels = "Total",
             indel_rate_sci = html("Rate (&times10<sup>&#8722<span>11</sup>)")) %>%
  cols_hide(columns = contains("E", ignore.case = F)) %>%
  tab_spanner(label = "Loss of Heterozygosity", 
              columns = matches("LOH|BY|BP|iLOH|tLOH")) %>%
  tab_spanner(label = "Point mutations", 
              columns = matches("SNM")) %>%
  tab_spanner(label = "Indels", 
              columns = matches("total_indels|indel_rate_sci")) %>%
  tab_footnote(footnote = html("per base-pair per generation &plusmn 
                               95% bootstrap confidence interval"),
               locations = cells_column_labels(columns = contains("rate"))) %>%
  # tab_footnote(footnote = html("&times 10<sup>-10</sup>"),
  #              locations = cells_column_labels(columns = c(BP_rate_sci, SNM_rate_sci))) %>%
  # tab_footnote(footnote = html("&times 10<sup>-11</sup>"),
  #              locations = cells_column_labels(columns = indel_rate_sci)) %>%
  tab_footnote(footnote = html("&plusmn binomial
                               95% confidence interval"),
               locations = cells_column_labels(columns = contains("f_"))) %>%
  cols_align(align = "right", columns = c(BY_label, iLOH_label)) %>%
  cols_width(Strain ~ px(80),
             contains("clones") ~ px(90),
             # matches("total&!BY") ~ px(80),
             # putative_LOH ~ px(90),
             # putative_BY ~ px(100),
             total_LOH ~ px(90),
             BY_label ~ px(140),
             # total_BY ~ px(90),
             contains("rate") ~ px(130),
             iLOH_label ~ px(120),
             # total_iLOH ~ px(90),
             # total_tLOH ~ px(90),
             # f_iLOH ~ px(130),
             total_SNM ~ px(60),
             total_indels ~ px(60)) %>%
  tab_style(style = cell_text(size = "large", weight = "bold"),
            locations =  cells_column_spanners())

master_table

# sum_px <- 3*80 + 2*90 + 5*100 + 2*120 + 150

gtsave(master_table, paste0(outIntDir, "master_table_2022_05.png"), zoom = 2,
       vwidth = 1800, vheight = 900)

# Expected and observed chromosomal LOH rates

# n clones
# chrom length
# n LOH expected
# n LOH observed

dip_length <- 2 * g_length
total_LOH <- sum(tx_LOHrate$total_LOH)
LOH_rate <- sum(tx_LOHrate$total_LOH)/sum(tx_LOHrate$n_clones)/n_gens/(2 * g_length)
LOH_rate_sci <- get_sci_digits(LOH_rate, base = -10)
LOH_boot <- boot(all_LOHcounts_merge_NS$n_LOH, mean.fun, R = 10000)

LOH_boot_CI <- boot.ci(boot.out = LOH_boot,
                      type = "basic")$basic[4:5]
LOH_rate_CI <- max(c(mean(all_LOHcounts_merge_NS$n_LOH) - LOH_boot_CI[1], 
      LOH_boot_CI[2] - mean(all_LOHcounts_merge_NS$n_LOH)))/n_gens/(2 * g_length)
LOH_CI_sci <- get_sci_digits(LOH_rate_CI, base = -10)

tx_LOHrate <- tx_LOHrate %>% mutate(bp_rate_sci = get_sci_digits(BP_rate))

expected_LOH_Tx <- n_clones_LOH_xTx %>%
  # ungroup() %>% 
  summarise(Tx_ID = Tx_ID, 
            expected_n = total_LOH * n / sum(n),
            expected_rate = total_LOH / sum(n) / n_gens / (2 * g_length),
            e_rate_sci = get_sci_digits(expected_rate)) %>%
  arrange(Tx_ID)

bp_CI_Tx <- tx_master_table %>% 
  select(BP_rate, BP_CIlo, BP_CIup) %>% 
  larger_CI() %>% get_sci_digits(base = -10)

# rom_CHROM chrom_length expected_LOH 
# total rate_sci rate_CI 
# expected_LOH_W total_LOH_W  e_bp_rate_W bp_rate_W_sci bp_CI_W
# expected_LOH_C total_LOH_C  e_bp_rate_C bp_rate_C_sci bp_CI_C
# expected_LOH_D total_LOH_D  e_bp_rate_D bp_rate_D_sci bp_CI_D 


tx_master_table %>% select(Strain, total_LOH, rate_sci, rate_CI)

wh_gen_row <- list("Genome", g_length, NA, 
                total_LOH, LOH_rate_sci, LOH_CI_sci,
                expected_LOH_Tx[1, 2], tx_LOHrate$total_LOH[1], expected_LOH_Tx[1, 4], tx_LOHrate$bp_rate_sci[1], bp_CI_Tx[1],
                expected_LOH_Tx[2, 2], tx_LOHrate$total_LOH[2], expected_LOH_Tx[2, 4], tx_LOHrate$bp_rate_sci[2], bp_CI_Tx[2],
                expected_LOH_Tx[3, 2], tx_LOHrate$total_LOH[3], expected_LOH_Tx[3, 4], tx_LOHrate$bp_rate_sci[3], bp_CI_Tx[3])


expected_TxCHROM_rate_wide <- expected_TxCHROM_rate %>% 
  select(Tx_ID, CHROM, chrom_length, expected_LOH, total_LOH, e_bp_rate, bp_rate, bp_CI_lo, bp_CI_up) %>%
  pivot_wider(id_cols = c(CHROM, chrom_length), names_from = Tx_ID, values_from = expected_LOH:bp_CI_up)
expected_TxCHROM_rate_wide$rom_CHROM <- levels(chrom_bound_BY$rom_CHROM)[
  as.numeric(expected_TxCHROM_rate_wide$CHROM)]


expected_TxCHROM_rate_wide <- expected_TxCHROM_rate_wide %>% 
  mutate(bp_rate_D_sci = get_sci_digits(bp_rate_D, base = -10),
         bp_rate_C_sci = get_sci_digits(bp_rate_C, base = -10),
         bp_rate_W_sci = get_sci_digits(bp_rate_W, base = -10),
         bp_rate_D_E = -10,
         bp_rate_C_E = -10,
         bp_rate_W_E = -10)

expected_TxCHROM_rate_wide$bp_CI_D <- round(larger_CI(expected_TxCHROM_rate_wide %>% 
                                               select(bp_rate_D, bp_CI_lo_D, bp_CI_up_D)) *
                                     10^-expected_TxCHROM_rate_wide$bp_rate_D_E, 2)
expected_TxCHROM_rate_wide$bp_CI_C <- round(larger_CI(expected_TxCHROM_rate_wide %>% 
                                                select(bp_rate_C, bp_CI_lo_C, bp_CI_up_C)) *
                                      10^-expected_TxCHROM_rate_wide$bp_rate_C_E, 2)
expected_TxCHROM_rate_wide$bp_CI_W <- round(larger_CI(expected_TxCHROM_rate_wide %>% 
                                                  select(bp_rate_W, bp_CI_lo_W, bp_CI_up_W)) *
                                        10^-expected_TxCHROM_rate_wide$bp_rate_W_E, 2)


LOHcounts_CHROM$rate_CI <- round(larger_CI(LOHcounts_CHROM %>% 
                                                        select(rate, rate_CI_95lo, rate_CI_95up)) *
                                              10^10, 2)
LOHcounts_CHROM$rate_sci <- get_sci_digits(LOHcounts_CHROM$rate, base = -10)


# expected_TxCHROM_rate_wide %>% select(contains("E", ignore.case = F))

expected_TxCHROM_rate_wide <- cbind(expected_TxCHROM_rate_wide, 
                                    LOHcounts_CHROM[, c("expected_LOH", "total", "rate_sci", "rate_CI")])


expected_TxCHROM_table_in <- 
  expected_TxCHROM_rate_wide %>% 
  select(rom_CHROM,
         chrom_length,
         expected_LOH,
         total,
         rate_sci,
         rate_CI,
         expected_LOH_W,
         total_LOH_W,
         e_bp_rate_W,
         bp_rate_W_sci,
         bp_CI_W,
         expected_LOH_C,
         total_LOH_C,
         e_bp_rate_C,
         bp_rate_C_sci,
         bp_CI_C,
         expected_LOH_D,
         total_LOH_D,
         e_bp_rate_D,
         bp_rate_D_sci,
         bp_CI_D)

expected_TxCHROM_table_in <- rbind(expected_TxCHROM_table_in, wh_gen_row)

expected_TxCHROM_table <- expected_TxCHROM_table_in %>% 
  # select(!contains("CIlo") & !contains("CIup")) %>%
  gt(rowname_col = "CHROM")  %>% 
  cols_merge_uncert(
    col_val = rate_sci,
    col_uncert = rate_CI) %>%
  cols_merge_uncert(
    col_val = bp_rate_W_sci,
    col_uncert = bp_CI_W) %>%
  cols_merge_uncert(
    col_val = bp_rate_C_sci,
    col_uncert = bp_CI_C) %>%
  cols_merge_uncert(
    col_val = bp_rate_D_sci,
    col_uncert = bp_CI_D) %>%
  fmt_number(columns = matches("expected_"), decimals = 1) %>%
  fmt_number(columns = matches("bp_"), decimals = 2) %>%
  cols_label(rom_CHROM = "Chr",
             chrom_length = "Chr length (bp)",
             expected_LOH = "Expected LOH",
             total = "Total LOH",
             rate_sci = html("LOH event rate (&times10<sup>&#8722<span>10</sup>)"),
             expected_LOH_W = "Expected LOH",
             total_LOH_W = "Total LOH",
             bp_rate_W_sci = html("LOH event rate (&times10<sup>&#8722<span>10</sup>)"),
             expected_LOH_C = "Expected LOH",
             total_LOH_C = "Total LOH",
             bp_rate_C_sci = html("LOH event rate (&times10<sup>&#8722<span>10</sup>)"),
             expected_LOH_D = "Expected LOH",
             total_LOH_D = "Total LOH",
             bp_rate_D_sci = html("LOH event rate (&times10<sup>&#8722<span>10</sup>)")) %>%
  cols_hide(columns = contains("e_bp", ignore.case = F)) %>%
  tab_spanner(label = "Total", 
              columns = c(expected_LOH, total, rate_sci)) %>%
  tab_spanner(label = "W Strain", 
              columns = c(expected_LOH_W, total_LOH_W, bp_rate_W_sci)) %>%
  tab_spanner(label = "C Strain", 
              columns = c(expected_LOH_C, total_LOH_C, bp_rate_C_sci)) %>%
  tab_spanner(label = "D Strain", 
              columns = c(expected_LOH_D, total_LOH_D, bp_rate_D_sci)) %>%
  tab_footnote(footnote = html("per base-pair per generation &plusmn 
                               95% bootstrap confidence interval, accounting for diploid genome"),
               locations = cells_column_labels(columns = contains("rate"))) %>%
  cols_width(rom_CHROM ~ px(80),
             chrom_length ~ px(100),
             contains("expected") ~ px(100),
             contains("total") ~ px(80),
             contains("rate") ~ px(120)) %>%
  tab_style(style = cell_text(size = "large", weight = "bold"),
            locations =  cells_column_spanners())

expected_TxCHROM_table

gtsave(expected_TxCHROM_table, paste0(outIntDir, "expected_TxCHROM_table.png"), 
       # zoom = 2,
       vwidth = 1600, vheight = 900)


