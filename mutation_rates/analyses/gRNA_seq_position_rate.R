# Determine whether there is enrichement for the gRNA sequence in the genome

###############################################################################
# Required functions and parameters ####

idx <- 1:length(gRNA_seq_split)
idx_w <- 1.46 + -0.112 * idx + 2.34E-03 * idx^2
idx_w <- ifelse(idx_w >1, 1, idx_w) - 1
get_score <- function(seq_in, gRNA_seq = gRNA_seq_split, 
                      idx = 1:20, idx_weights = idx_w, run_weight = 1) {
  # linear model with position weighting 0.125 - 1 5'-3'
  # and weighted mismatch runs a(sum(i:j)) + b(length)
  # seq_in <- gRNA_seq
  # seq_in <- chr_seq_v[27224:27243]
  # seq_in <- c(gRNA_seq[11:20], gRNA_seq[11:20])
  # seq_in <- c(gRNA_seq[1:6], gRNA_seq[3:7], gRNA_seq[12:20])
  # p = 1.5061 * exp(-0.1137 * i)
  # relative_activity = 1.46 + -0.112x + 2.34E-03x^2
  # max single mismatch = 25% activity
  # abs_activity = 0.5/1.288 * total_penalty + 1
  i_match <- gRNA_seq == seq_in
  rle_match <- rle(i_match) %>% unclass() %>% as.data.frame()
  i_rle_mis <- !rle_match$values
  mid_i_match <- cumsum(rle_match$lengths) - (rle_match$lengths - 1)/2
  i_weights <- matrix(c(ifelse(floor(mid_i_match) > 0, floor(mid_i_match), 1), 
                        ifelse(ceiling(mid_i_match) > 0, ceiling(mid_i_match), 1)), 
                        ncol = 2)
  penalty_df <- data.frame(match = !i_rle_mis, length = rle_match$lengths,
                           weight = apply(i_weights, MARGIN = 1, function(x) mean(idx_weights[x])))
  penalty_df$total <- penalty_df$length * penalty_df$weight
  final_score <- sum(penalty_df$total[i_rle_mis] + 
    run_weight * (penalty_df$length[i_rle_mis] - 1) * penalty_df$weight[i_rle_mis])
  
  return(final_score)
}


###############################################################################
# Prepare reference and calculate nucleotide frequencies ####
# seq_in <- BY4742_seq
# seq_v_in <- BY4742_seq_v
# seq_v_in <- RM_seq_v

nuc_counts <- lapply(seq_v_in, function(x) seqinr::count(x, wordsize = 1))
nuc_df <- data.frame(a = 0, c = 0, g = 0, t = 0)
for(chr in 1:(length(nuc_counts) - 1)) {
  nuc_df <- nuc_counts[[chr]] + nuc_df
}

nuc_freq_df <- t(nuc_df) %>% as.data.frame() %>% rownames_to_column()
names(nuc_freq_df) <- c("nuc", "count")
nuc_freq_df <- nuc_freq_df %>% mutate(freq = count/sum(count))

###############################################################################
# Calculate gRNA similarity scores for reference sequence ####

pams <- list("ngg" = "gg", "nag" = "ag", "ngg_rev" = "cc", "nag_rev" = "ct")
genome_scores <- data.frame(NULL)
for(chr in 1:16) {
  # chr <- 12
  chr_seq <- as.character(seq_in[chr])
  # str_locate(chr_seq, gRNA_seq)
  chr_seq_v <- seq_v_in[[chr]]
  chrom_pams <- lapply(pams, function(x) str_locate_all(chr_seq, x)[[1]][, 2])
  chrom_pams <- lapply(chrom_pams, function(x) x[x > 22 & x < (chrom_lengths_BY[chr] - 22)])
  chrom_scores <- list()
  for(p in 1:2) {
    # p <- 1
    # str_locate_all returns the index of the end of the PAM, 
    # indicies for target sequence are 3bp and 22bp upstream
    align_scores <- sapply(chrom_pams[[p]], function(x) get_score(chr_seq_v[(x - 22):(x - 3)]))
    chrom_scores[[p]] <- align_scores
  }
  for(p in 3:4) {
    # p <- 3
    # str_locate_all returns the index of the end of the PAM in the 
    # forward direction, which is the first position of the rev comp PAM
    # indicies for rev comp target sequence are 1bp and 20bp downstream
    # also, the seq and gRNA must be reversed to index correctly in get_score()
    align_scores <- sapply(chrom_pams[[p]], 
                           function(x) get_score(chr_seq_v[(x + 2):(x + 21)],
                                                 gRNA_seq = rev(comp(gRNA_seq_comp_split))))
    chrom_scores[[p]] <- align_scores
  }
  names(chrom_scores) <- names(pams)
  chrom_scores_df <- stack(chrom_scores)
  colnames(chrom_scores_df) <- c("score", "pam")
  chrom_scores_df$POS <- unlist(chrom_pams,  use.names = FALSE)
  chrom_scores_df$CHROM <- chrom_bound_BY$CHROM[chr]
  genome_scores <- rbind(genome_scores, chrom_scores_df)
}

genome_scores <- data.frame(NULL)
  for(chr in 1:16) {
    # chr <- 4
    chr_seq <- as.character(BY4742_seq[chr])
    # str_locate(chr_seq, gRNA_seq)
    chr_seq_v <- BY4742_seq_v[[chr]]
    chrom_pams <- lapply(pams, function(x) str_locate_all(chr_seq, x)[[1]][, 2])
    chrom_pams <- lapply(chrom_pams, function(x) x[x > 27 & x < (chrom_lengths_BY[chr] - 27)])
    chrom_scores <- list()
    for(p in 1:2) {
      # p <- 1
      # str_locate_all returns the index of the end of the PAM, 
      # indicies for target sequence are 3bp and 22bp upstream
      align_scores <- lapply(chrom_pams[[p]], 
                             function(x) get_score_NW(patt = gRNA_seq_split,
                                                      subj = chr_seq_v[(x - 27):(x - 3)],
                                                      params = NW_params)) %>% 
        bind_rows(.id = "pam") %>% as.data.frame()
      chrom_scores[[p]] <- align_scores
    }
    for(p in 3:4) {
      # p <- 3
      # str_locate_all returns the index of the end of the PAM in the 
      # forward direction, which is the first position of the rev comp PAM
      # indicies for rev comp target sequence are 1bp and 20bp downstream
      # also, the seq and gRNA must be reversed to index correctly in get_score()
      align_scores <- lapply(chrom_pams[[p]], 
                             function(x) get_score_NW(patt = rev(comp(gRNA_seq_split)),
                                                      subj = chr_seq_v[(x + 2):(x + 26)],
                                                      reverse = F,
                                                      params = NW_params)) %>% 
        bind_rows(.id = "pam") %>% as.data.frame()
      chrom_scores[[p]] <- align_scores
    }
    names(chrom_scores) <- names(pams)
    chrom_scores_df <- bind_rows(chrom_scores, .id = "pam")
    # colnames(chrom_scores_df) <- c("score", "pam")
    chrom_scores_df$POS <- unlist(chrom_pams,  use.names = FALSE)
    chrom_scores_df$CHROM <- chrom_bound_BY$CHROM[chr]
    chrom_scores_df <- chrom_scores_df %>% filter(score >= -5)
    genome_scores <- rbind(genome_scores, chrom_scores_df)
}

chr_seq_v[27224:27243]
gRNA_seq_comp_split
chr_seq_v[1028670:1028689]
gRNA_seq_split

# genome_scores_BY <- genome_scores
# genome_scores_RM <- genome_scores
genome_scores %>% ggplot() + geom_histogram(aes(x = score), binwidth = 0.1)

genome_scores %>% filter(score > -1) %>% ggplot() + geom_freqpoly(aes(x = score), binwidth = 0.05)
genome_scores %>% filter(score > -1) %>% arrange(desc(score))

genome_scores_RM$POSi <- 0
genome_scores_RM_lift <- RMxBY_liftover(genome_scores_RM)

genome_scores_all <- merge(genome_scores_BY, 
                           genome_scores_RM_lift, 
                           by = c("CHROM", "POS", "pam"), 
                           suffixes = c("_BY", "_RM"), 
                           all = T)

genome_scores_all$score_RM[is.na(genome_scores_all$score_RM)] <-
  genome_scores_all$score_BY[is.na(genome_scores_all$score_RM)]
genome_scores_all$score_BY[is.na(genome_scores_all$score_BY)] <-
  genome_scores_all$score_RM[is.na(genome_scores_all$score_BY)]

genome_scores_all$score_max <- apply(genome_scores_all[, c("score_BY", "score_RM")], 
                                     MARGIN = 1, 
                                     function(x) max(x))

# cs <- colsplit(genome_scores$pam, pattern = "_", names = c("pam", "dir"))
# genome_scores <- genome_scores %>% select(-pam) %>% mutate(pam = cs$pam, dir = cs$dir)
# genome_scores$dir[genome_scores$dir == ""] <- "fwd"
# genome_scores$exp_activity <- ifelse(1 + genome_scores$score < 0, 0, 1 + genome_scores$score)
# genome_scores$exp_activity <- ifelse(genome_scores$pam == "ngg", genome_scores$exp_activity, genome_scores$exp_activity * 0.5)

genome_scores_all %>% ggplot() + geom_histogram(aes(x = score_BY))
genome_scores_all %>% ggplot() + geom_histogram(aes(x = score_RM))
genome_scores_all %>% ggplot() + geom_histogram(aes(x = score_max))

###############################################################################
# Test whether the similarity score is higher at D hotspots ####

i_hs1 <- genome_scores_all$CHROM == "05" & 
  genome_scores_all$POS >= 180000 & 
  genome_scores_all$POS <= 271000
genome_scores_all$hs1 <- ifelse(i_hs1, "y", "n")
genome_scores_all$hs1 <- factor(genome_scores_all$hs1)

i_hs2 <- genome_scores_all$CHROM == "02" & 
  genome_scores_all$POS >= 325000 & 
  genome_scores_all$POS <= 365000
genome_scores_all$hs2 <- ifelse(i_hs2, "y", "n")
genome_scores_all$hs2 <- factor(genome_scores_all$hs2)

group_SW_hotspots <- function(x, gap = 100000) {
  # x <- PM_BPrate_SW %>% filter(pm > 5, Tx_ID == "D") %>% pull(start)
  d_x <- x - c(x[1], x[1:(length(x) - 1)])
  i_new <- c(0, which(d_x > gap), length(x))
  g <- 1
  group <- rep(0, length(x))
  for(i in 2:length(i_new)) {
    # i <- 8
    group[(i_new[i - 1]):i_new[i]] <- g
    g <- g + 1
  }
  return(group)
}

pm_hs_POSi <- PM_BPrate_SW %>% 
  filter(pm > 5, Tx_ID == "D") %>% 
  mutate(hs = group_SW_hotspots(start)) %>%
  group_by(hs) %>%
  summarize(hs_start = min(start), hs_end = max(end))

genome_scores_all$hs_PM <- "n"
for(i in 1:nrow(pm_hs_POSi)) {
  i_hs <- genome_scores_all$POSi > pm_hs_POSi$hs_start[i] &
    genome_scores_all$POSi < pm_hs_POSi$hs_end[i]
  genome_scores_all$hs_PM[i_hs] <- "y"
}

genome_scores_all %>% ggplot() + geom_histogram(aes(x = score_max)) + facet_grid(hs2~., scales = "free_y")
genome_scores_all %>% filter(score_max < 0) %>% group_by(hs2) %>% summarize(m = max(score_max, na.rm = T))

hs1_perm_mean <- perm_test(genome_scores_all, 
                           cat_var = "hs1", cat_names = c("n", "y"), 
                           response_var = "score_max", 
                           test_stat = mean)

hs1_perm_max <- perm_test(genome_scores_all, 
                          cat_var = "hs1", cat_names = c("n", "y"), 
                          response_var = "score_max", 
                          test_stat = max)

hs2_perm_mean <- perm_test(genome_scores_all, 
                           cat_var = "hs2", cat_names = c("n", "y"), 
                           response_var = "score_max", 
                           test_stat = mean)

hs2_perm_max <- perm_test(genome_scores_all, 
                          cat_var = "hs2", cat_names = c("n", "y"), 
                          response_var = "score_max", 
                          test_stat = max)

hs_PM_perm_mean <- perm_test(genome_scores_all, 
                           cat_var = "hs_PM", cat_names = c("n", "y"), 
                           response_var = "score_max", 
                           test_stat = mean, n_perms = 100)

hs_PM_perm_max <- perm_test(genome_scores_all, 
                          cat_var = "hs_PM", cat_names = c("n", "y"), 
                          response_var = "score_max", 
                          test_stat = max, n_perms = 100)

###############################################################################
# Test whether there is a higher density of PAMs at C hotspots ####
min_score_g <- -10
genome_scores_all$CHROM <- factor(genome_scores$CHROM)
genome_scores_all$POSi <- genome_scores %>% ConvertPosIndicies(pos_col = "POS", index_out = "POSi")

pam_density_50kb <- genome_scores_all %>% 
  filter(score_max >= min_score_g) %>%
  mutate(p = 1) %>%
  SliderCalc(data_col = "p", index_col = "POSi", window_size = 50000, 
             slide_interval = 5000, chrom_win = T, summary_stat = sum)

# genome_scores_50kb$start_POS <- genome_scores_50kb %>% ConvertPosIndicies(pos_col = "start", index_out = "POS", add_chroms = F)
# genome_scores_50kb <- cbind(genome_scores_50kb, start_POS = genome_scores_POS$POS, CHROM = genome_scores_POS$CHROM)
# genome_scores_50kb$end_POS <- genome_scores_50kb %>% ConvertPosIndicies(pos_col = "end", index_out = "POS")
# genome_scores_50kb <- genome_scores_50kb %>% mutate(mid_POS = round((start_POS + end_POS)/2))

LOH_BPrate_SW_C <- LOH_BPrate_SW %>% 
  filter(Tx_ID == "C")
LOH_BPrate_SW_C <- PM_BPrate_SW %>% 
  mutate(bp = pm) %>%
  # filter(pm > 5) %>%
  filter(Tx_ID == "C")

LOH_BPrate_SW_C <- merge(LOH_BPrate_SW_C, pam_density_50kb[, c("start", "p")], by = "start", all = T)

LOH_BPrate_SW_C %>% 
  # filter(CHROM == "07") %>%
  ggplot() + geom_point(aes(x = bp, y = p)) + facet_wrap(~CHROM, ncol = 4)

LOH_BPrate_SW_C %>% 
  mutate(bp = pm) %>%
  # filter(CHROM != "09") %>% 
  ggplot() + geom_histogram(aes(x = bp))

LOH_BPrate_SW_C %>% filter(CHROM != "09") %>% filter(bp > 5)

i_C_hs <- LOH_BPrate_SW_C$bp > 5
i_C_hs <- LOH_BPrate_SW_C$pm > 5
LOH_BPrate_SW_C$hs_C <- ifelse(i_C_hs, "y", "n")

# i_C_hs <- LOH_BPrate_SW_C$CHROM != "09" & LOH_BPrate_SW_C$bp > 5
# LOH_BPrate_SW_C$hs_C <- ifelse(i_C_hs, "y", "n")

hs2_perm_mean <- perm_test(LOH_BPrate_SW_C, 
                           cat_var = "hs_C", cat_names = c("n", "y"), 
                           response_var = "p", 
                           test_stat = mean)

hs2_perm_max <- perm_test(LOH_BPrate_SW_C, 
                           cat_var = "hs_C", cat_names = c("n", "y"), 
                           response_var = "p", 
                           test_stat = max)

###############################################################################
# Distribution of scores for random sequences ####
nucs <- c("a", "c", "g", "t")
n_seqs <- nrow(genome_scores)
set.seed(313)
seqs <- matrix(sample(nucs, 20 * n_seqs, replace = T, prob = nuc_freq_df$freq), ncol = 20)
seqs <- rbind(gRNA_seq_split, seqs)
# for each alignment, get indicies of matches
# random_indices <- apply(seqs, MARGIN = 1, function(x) c(1:20)[gRNA_seq_split == x])
# random_scores <- lapply(random_indices, function(x) get_score(x))
random_scores <- apply(seqs, MARGIN = 1, function(x) get_score(x, run_weight = 1))

# Check top scores in genome
top_hits <- genome_scores %>% 
  arrange(desc(score)) %>%
  dplyr::slice(1:10)

# Since there are no other perfect matches aside from the ade2 target, 
# label target by score of 0
genomeVsRandom_scores$target <- F
genomeVsRandom_scores$target[genomeVsRandom_scores$score == 0] <- T

genomeVsRandom_scores <- genome_scores %>% filter(!target) %>% arrange(score)
genomeVsRandom_scores$random <- sort(random_scores[2:(nrow(genome_scores) + 1)])
genomeVsRandom_scores$pam <- factor(genomeVsRandom_scores$pam, levels = c("ngg", "nag"))
genomeVsRandom_scores$PAM <- factor(genomeVsRandom_scores$pam, labels = c("NGG", "NAG"))


# QQplot of observed sequence similarity scores vs scores for random sequences
###############################################################################

genomeVsRandom_plot <- genomeVsRandom_scores %>% 
  filter(!target) %>%
  sample_n(10000, replace = F) %>%
  ggplot() + 
  geom_point(aes(x = random, y = score, color = PAM), size = 2) +
  geom_abline(color = "blue3") +
  xlab("Random sequences score") + 
  ylab("PAM-adjacent sequence score") +
  # facet_grid(~PAM)
  scale_color_manual(values = c("grey50", "black")) +
  theme(legend.position = c(0.9, 0.1),
        legend.background = element_rect(fill = "white", color = "white"),
        text = element_text(size = 18))

genomeVsRandom_plot

ggsave(file.path(outIntDir, "genomeVsRandom_plot.png"), 
       plot = genomeVsRandom_plot,
       device = "png",
       width = 11, height = 8.5, 
       units = "in",
       dpi = 600)

df <- data.frame(i = 0:n_seqs, score = unlist(random_scores))
# df_g <- data.frame(i = 0, score = get_score(gRNA_seq_split))
# df <- rbind(df_g, df)
genome_scores %>% 
  chrom_scores_df %>%
  # df %>%
  # filter(score >= 0) %>%
  ggplot + geom_histogram(aes(x = score)) +
  xlim(NA, 0)

min_score_g <- genome_scores %>% arrange(desc(score)) %>%
  dplyr::slice(1:ceiling(0.01 * nrow(genome_scores))) %>% 
  pull(score) %>% min()

genome_scores %>% 
  filter(score >= min_score_g) %>%
  # arrange(score)
  ggplot() + 
  # geom_line(aes(x = POS, y = score)) +
  geom_point(aes(x = POS, y = score, color = pam)) +
  facet_wrap(~CHROM, ncol = 4, scales = "free_x")


genome_scores$CHROM <- factor(genome_scores$CHROM)
genome_scores$POSi <- genome_scores %>% ConvertPosIndicies(pos_col = "POS", index_out = "POSi")

genome_scores_50kb <- genome_scores %>% 
  filter(score >= min_score_g) %>%
  SliderCalc(data_col = "score", index_col = "POSi", window_size = 50000, 
             slide_interval = 50000, chrom_win = T, summary_stat = max)

genome_scores_50kb$start_POS <- genome_scores_50kb %>% ConvertPosIndicies(pos_col = "start", index_out = "POS", add_chroms = F)
# genome_scores_50kb <- cbind(genome_scores_50kb, start_POS = genome_scores_POS$POS, CHROM = genome_scores_POS$CHROM)
genome_scores_50kb$end_POS <- genome_scores_50kb %>% ConvertPosIndicies(pos_col = "end", index_out = "POS")

genome_scores_50kb <- genome_scores_50kb %>% mutate(mid_POS = round((start_POS + end_POS)/2))

LOH_rate_SW_wide <- LOH_mid_SW %>% 
  filter(Tx_name != "Cas9") %>% 
  select(c(Tx_name, mid_POSi, rate)) %>% 
  pivot_wider(names_from = Tx_name, values_from = rate)

LOH_rate_SW_wide <- LOH_rate_SW_wide %>% mutate(rel_rate = (Drive - WT)/(Drive + WT))

LOH_rate_SW_wide$score <- genome_scores_50kb$score
LOH_rate_SW_wide$target <- F
LOH_rate_SW_wide$target[LOH_rate_SW_wide$score == 0] <- T

LOH_rate_SW_wide %>% 
  filter(WT > 0 | Drive > 0) %>%
  ggplot() + geom_abline(color = "grey30") + 
  geom_jitter(aes(x = WT, y = Drive, color = score)) +
  scale_color_gradient2(low = "grey70", mid = "skyblue", high = "red3", midpoint = -1.5)

score_lm <- lm(rel_rate ~ score, data = LOH_rate_SW_wide %>% filter(target == F))
summary(score_lm)
lm_b <- score_lm$coefficients[1]
lm_m <- score_lm$coefficients[2]

segment_coord <- c(x = min(LOH_rate_SW_wide$score), 
                   xend = max(LOH_rate_SW_wide$score), 
                   y = min(LOH_rate_SW_wide$score) * lm_m + lm_b, 
                   yend = max(LOH_rate_SW_wide$score) * lm_m + lm_b)

# GET EQUATION AND R-SQUARED AS STRING
# SOURCE: https://groups.google.com/forum/#!topic/ggplot2/1TgH-kG5XMA

lm_eqn <- function(lm_in){
  eq <- substitute(italic(y) == a + b * italic(x)*","~~italic(r)^2~"="~r2, 
                   list(a = format(unname(coef(lm_in)[1]), digits = 2),
                        b = format(unname(coef(lm_in)[2]), digits = 2),
                        r2 = format(summary(lm_in)$r.squared, digits = 3, scientific = T)))
  as.character(as.expression(eq))
}


# Plot relationship of similarity score and rate of LOH in drive above WT #####
###############################################################################
LOHrate_vs_seqScore <- LOH_rate_SW_wide %>%
  filter(!is.na(rel_rate), !target) %>%
  ggplot() +
  geom_point(aes(x = score, y = rel_rate, color = target)) +
  annotate(geom = "segment", x = segment_coord[1], xend = segment_coord[2], 
           y = segment_coord[3], yend = segment_coord[4]) + 
  annotate(geom = "label", x = -0.5, y = 0.75, label = lm_eqn(score_lm), parse = TRUE, 
           label.size = 0, label.padding = unit(0.5, "lines")) +
  scale_color_manual(values = c("black", "red3"), guide = NULL) +
  # geom_text(x = -0.75, y = 0.75, label = lm_eqn(score_lm), parse = TRUE) +
  xlab("gRNA sequence similarity score") +
  ylab("Normalized LOH rate difference in Drive vs WT") +
  xlim(NA, 0) +
  theme(text = element_text(size = 16))

LOHrate_vs_seqScore

ggsave(file.path(outIntDir, "LOHrate_vs_seqScore_plot.png"), 
       plot = LOHrate_vs_seqScore,
       device = "png",
       width = 11, height = 8.5, 
       units = "in",
       dpi = 600)


# distribution of scores for random sequences
nucs <- c("a", "t", "c", "g")
n_seqs <- 100
set.seed(313)
seqs <- matrix(sample(nucs, 20 * n_seqs, replace = T), ncol = 20)

random_scores <- lapply(seqs, function(x) get_score(x))

df <- data.frame(i = 1:n_seqs, score = unlist(random_scores))

min_score <- df %>% arrange(desc(score)) %>%
  dplyr::slice(1:ceiling(0.01 * nrow(df))) %>% 
  pull(score) %>% min()

df %>% ggplot() + geom_histogram(aes(x = score)) +
  # xlim(0, 210) +
  geom_vline(aes(xintercept = min_score), color = "blue2")

df %>% ggplot() + geom_histogram(aes(x = log(score, 10)), binwidth = 0.01) # +
  geom_vline(aes(xintercept = log(min_score, 10)), color = "blue2")

genome_scores_all %>% ggplot() + geom_histogram(aes(x = score), binwidth = 1)
  
genome_scores_all %>% ggplot() + geom_histogram(aes(x = log(score, 10)), binwidth = 0.01) +
  geom_vline(aes(xintercept = log(min_score, 10)), color = "blue2")

genome_scores %>%
  filter(CHROM == "05") %>% 
  # filter(score > 300) %>% 
  ggplot() + 
  # geom_line(aes(x = POS, y = log(score, 10)))
  geom_line(aes(x = POS, y = score))

# genome_scores_rev %>% arrange(desc(score)) %>% dplyr::slice(1:ceiling(0.01 * nrow(df))) %>% pull(score) %>% min()

genome_scores %>% 
  # filter(score >= min_score) %>%
  ggplot() + geom_histogram(aes(x = log(score, 10)), binwidth = 0.01) +
  geom_vline(aes(xintercept = log(min_score, 10)), color = "blue2")

genome_scores_all %>% 
  filter(score >= 1) %>%
  ggplot() + geom_histogram(aes(x = log(score, 10)))

chr <- "05"
genome_scores_all %>% 
  filter(score >= min_score) %>%
  # filter(CHROM == chr) %>%
  ggplot() + geom_line(aes(x = POS, y = score)) + 
  # geom_point(data = gRNA_hits %>% filter(hit_score > 5, CHROM == chr), 
  #           aes(x = s_start, y = hit_score), color = "blue2") +
  facet_wrap(~CHROM, ncol = 4, scales = "free_x") +
  ylim(min_score, 210)

genome_scores_all %>% 
  filter(score >= 140) %>%
  filter(CHROM == chr)

genome_scores_all %>% 
  arrange(desc(score)) %>%
  dplyr::slice(1:10)

c <- 15
POS <- 565935
chr_seq_v <- BY_seq_v[[c]]

chr_seq_v[POS:(POS + 22)]
gRNA_seq_split

chr_seq_v[(POS - 22):POS]
gRNA_seq_comp_split


genome_scores_all %>% 
  filter(score >= min_score) %>%
  filter(CHROM == "15", POS > 550000, POS < 600000) %>%
  summarize(av = mean(score))

genome_scores_all %>% 
  filter(score >= min_score) %>%
  filter(CHROM == "15", POS > 500000, POS < 550000) %>%
  summarize(av = mean(score))




# Density of PAM sites among chromosomes ######################################

n_PAMs_ngg <- data.frame(CHROM = chrom_bound_BY$CHROM,  n_pam = 0)
genome_pams_ngg <- data.frame(NULL)
for(c in 1:16) {
  # c <- 7
  chr_seq <- as.character(BY_seq[c])
  chr_seq_v <- BY_seq_v[[c]]
  
  chrom_pams_ngg <- str_locate_all(chr_seq, "gg")[[1]][, 2]
  chrom_pams_ccn <- str_locate_all(chr_seq, "cc")[[1]][, 2]
  chrom_pams_ngg <- chrom_pams_ngg[chrom_pams_ngg > 22]
  chrom_pams_ccn <- chrom_pams_ccn[chrom_pams_ccn > 22]
  # n_PAMs_ngg$n_pam[c] <- length(chrom_pams_ngg) + length(chrom_pams_ccn)
  chrom_pams <- data.frame(CHROM = chrom_bound_BY$CHROM[c], POS = c(chrom_pams_ngg, chrom_pams_ccn))
  genome_pams_ngg <- rbind(genome_pams_ngg, chrom_pams)
}
n_PAMs_ngg$rate <- n_PAMs_ngg$n_pam/chrom_lengths_BY
n_PAMs_ngg$seq <- "NGG"

n_PAMs_nga <- data.frame(CHROM = chrom_bound_BY$CHROM,  n_pam = 0)
genome_pams_nga <- data.frame(NULL)
for(c in 1:16) {
  # c <- 16
  chr_seq <- as.character(BY_seq[c])
  chr_seq_v <- BY_seq_v[[c]]
  
  chrom_pams_nga <- str_locate_all(chr_seq, "ga")[[1]][, 2]
  chrom_pams_tcn <- str_locate_all(chr_seq, "tc")[[1]][, 2]
  chrom_pams_nga <- chrom_pams_nga[chrom_pams_nga > 22]
  chrom_pams_tcn <- chrom_pams_tcn[chrom_pams_tcn > 22]
  # n_PAMs_nga$n_pam[c] <- length(chrom_pams_nga) + length(chrom_pams_tcn)
  chrom_pams <- data.frame(CHROM = chrom_bound_BY$CHROM[c], POS = c(chrom_pams_nga, chrom_pams_tcn))
  genome_pams_nga <- rbind(genome_pams_nga, chrom_pams)
}

n_PAMs_nga$rate <- n_PAMs_nga$n_pam/chrom_lengths_BY
n_PAMs_nga$seq <- "NGA"

n_PAMs <- rbind(n_PAMs_ngg, n_PAMs_nga)
n_PAMs$seq <- factor(n_PAMs$seq)

n_PAMs %>% ggplot() + geom_col(aes(x = CHROM, y = rate, fill = seq))

genome_pams_ngg$seq <- "NGG"
genome_pams_ngg$POSi <- ConvertPosIndicies(genome_pams_ngg)
genome_pams_ngg$p <- 1

genome_pams_nga$seq <- "NGA"
genome_pams_nga$POSi <- ConvertPosIndicies(genome_pams_nga)
genome_pams_nga$p <- 1

genome_pams <- rbind(genome_pams_ngg, genome_pams_nga)
genome_pams$seq <- factor(genome_pams$seq)

genome_pams_50k <- SliderCalc(genome_pams, data_col = "p", index_col = "POSi", factor_col = "seq",
                                  window_size = 50000, chrom_win = T, summary_stat = sum)

genome_pams_ngg_50k <- SliderCalc(genome_pams_ngg, data_col = "p", index_col = "POSi",
                              window_size = 50000, chrom_win = T, summary_stat = sum)

genome_pams_nga_50k <- SliderCalc(genome_pams_nga, data_col = "p", index_col = "POSi",
                                  window_size = 50000, chrom_win = T, summary_stat = sum)

genome_pams_50k_plot <- genome_pams_50k %>% 
  ggplot() + 
  geom_line(aes(x = (start + end)/2, y = p, color = seq)) +
  facet_wrap(~CHROM, scales = "free_x", ncol = 4) +
  scale_color_manual(values = c("darkred", "darkblue")) +
  xlab("Window position") + ylab("Number of PAMs")

ggsave(file.path(outIntDir, "genome_pams_50k_plot.png"), 
       plot = genome_pams_50k_plot,
       device = "png",
       width = 16, height = 9, 
       units = "in",
       dpi = 600)


# BLAST gRNA method. NOT USED #################################################

score_gRNA_hits <- function(b_table, gRNA_len = 20) {
  # b_table <- blast_table
  end_score <- (b_table$q_end - gRNA_len)*4
  start_score <- 1 - b_table$q_start
  match_score <- 0 - b_table$mismatch*2
  gap_score <- 0 - b_table$gaps*3
  total_score <- gRNA_len + end_score + start_score + match_score + gap_score
}

importBlastTables <- function(file_list) {
  cn <- c("Job", "ID", "identity", "length", 
          "mismatch", "gaps", "q_start", "q_end",
          "s_start", "s_end", "e_value", "bit_score")
  table_list <- list(NULL)
  for(f in 1:length(file_list)) {
    # f <- 1
    fn <- blast_list[[f]]
    blast_table <- read_delim(fn, "\t", col_names = F, comment = "#") %>%
      as.data.frame()
    
    colnames(blast_table) <- cn
    
    blast_table <- blast_table %>% select(!Job)
    
    blast_table <- blast_table %>% filter(ID %in% chrom_IDs$ID)
    
    blast_table$CHROM <- unlist(lapply(blast_table$ID, 
                                       function(x) chrom_IDs$CHROM[chrom_IDs$ID == x]))
    
    blast_table$s_start_POSi <- ConvertPosIndicies(blast_table, pos_col = "s_start", index_out = "POSi")
    blast_table$s_end_POSi <- ConvertPosIndicies(blast_table, pos_col = "s_end", index_out = "POSi")
    
    blast_table$strand <- "plus"
    iMinus <- blast_table$s_end_POSi < blast_table$s_start_POSi
    blast_table$strand[iMinus] <- "minus"
    # bt_2 <- blast_table
    m_end <- blast_table$s_start[iMinus]
    blast_table$s_start[iMinus] <- blast_table$s_end[iMinus]
    blast_table$s_end[iMinus] <- m_end
    
    m_end <- blast_table$s_start_POSi[iMinus]
    blast_table$s_start_POSi[iMinus] <- blast_table$s_end_POSi[iMinus]
    blast_table$s_end_POSi[iMinus] <- m_end
    blast_table <- blast_table %>% arrange(s_start_POSi)
    
    blast_table$mid_POSi <- round(blast_table$s_start + blast_table$s_end)/2
    blast_table$hit_score <- score_gRNA_hits(blast_table)
    table_list[[f]] <- blast_table
  }
  names(table_list) <- paste0("table_", 1:length(file_list))
  return(table_list)
}

getUniqueHits <- function(blast_list) {
  # blast_list <- gRNA_blast_tables
  l <- length(blast_list)
  c_1 <- unlist(sapply(1:l, function(x) rep(x, l - x)))
  c_2 <- unlist(sapply(1:(l-1), function(x) seq.int(x+1, l)))
  l_combos <- data.frame(c_1, c_2)
  blast_table_1 <- blast_list[[1]]
  for(i in 2:l) {
    # i <- 3
    # blast_table_1 <- blast_list[[l_combos[i, 1]]]
    blast_table_2 <- blast_list[[i]]
    # blast_table_2_2 <- blast_table_2
    iDup_s_POSi <- blast_table_2$s_start_POSi %in% blast_table_1$s_start_POSi
    iDup_e_POSi <- blast_table_2$s_end_POSi %in% blast_table_1$s_end_POSi
    blast_table_2 <- blast_table_2[!(iDup_s_POSi & iDup_e_POSi), ]
    
    blast_ranges_1 <- IRanges(start = blast_table_1$s_start_POSi, end = blast_table_1$s_end_POSi)
    blast_ranges_2 <- IRanges(start = blast_table_2$s_start_POSi, end = blast_table_2$s_end_POSi)
    # overlap_1_1 <- findOverlaps(blast_ranges_1, drop.self = T, select = "first")
    # overlap_2_2 <- findOverlaps(blast_ranges_2, drop.self = T, select = "first")
    
    overlap_1_2 <- findOverlaps(blast_ranges_1, blast_ranges_2, select = "first")
    overlap_valid_1_2 <- overlap_1_2[!is.na(overlap_1_2)]
    
    Unique_1 <- blast_table_1[is.na(overlap_1_2), ]
    Overlap_1 <- blast_table_1[!is.na(overlap_1_2), ]
    Unique_2 <- blast_table_2[-overlap_valid_1_2, ]
    Overlap_2 <- blast_table_2[overlap_valid_1_2, ]
    
    scores_1 <- Overlap_1$hit_score
    scores_2 <- Overlap_2$hit_score
    q_ends_1 <- Overlap_1$q_end
    q_ends_2 <- Overlap_2$q_end
    i_1 <- (scores_1 > scores_2 | (scores_1 == scores_2 & q_ends_1 >= q_ends_2))
    i_2 <- (scores_2 > scores_1 | (scores_1 == scores_2 & q_ends_1 < q_ends_2))
    best_algn_1 <- Overlap_1[i_1, ]
    best_algn_2 <- blast_table_2[overlap_valid_1_2[i_2], ]
    combo_out <- rbind(Unique_1, Unique_2, best_algn_1, best_algn_2) %>% arrange(s_start_POSi)
    rownames(combo_out) <- NULL
    
    dup_s_POSi <- combo_out$s_start_POSi[base::duplicated(combo_out$s_start_POSi)]
    dup_e_POSi <- combo_out$s_end_POSi[base::duplicated(combo_out$s_end_POSi)]
    dups <- combo_out %>% filter(s_start_POSi %in% dup_s_POSi, s_end_POSi %in% dup_e_POSi) %>% nrow()
    # if(dups != 0) {
    #   print("duplicates found")
    # }
    combo_out <- combo_out %>% distinct(s_start_POSi, s_end_POSi, .keep_all = T)
    blast_table_1 <- combo_out
  }
  return(blast_table_1)
}

gRNA_blast_tables <- importBlastTables(blast_list)
gRNA_hits <- getUniqueHits(gRNA_blast_tables)

min(gRNA_hits$length)

gRNA_hits %>% ggplot() + geom_histogram(aes(x = length), binwidth = 1) + xlim(0, NA)
gRNA_hits %>% ggplot() + geom_histogram(aes(x = q_end), binwidth = 1) + xlim(0, NA)
gRNA_hits %>% ggplot() + geom_histogram(aes(x = q_start), binwidth = 1) + xlim(0, 20)

gRNA_hits %>% filter(q_end >= 17) %>% ggplot() + geom_bar(aes(x = CHROM))


gRNA_hits %>% ggplot() + geom_point(aes(x = log10(1/e_value), y = hit_score))
gRNA_hits %>% ggplot() + geom_jitter(aes(x = q_end, y = hit_score))

gRNA_hits %>% 
  filter(q_end >= 16) %>%
  ggplot() + geom_histogram(aes(x = hit_score), binwidth = 1) +
  facet_grid(q_end~.)

gRNA_hits %>% 
  filter(CHROM == "05") %>%
  # filter(q_end >= 16) %>%
  filter(hit_score > 0) %>% 
  ggplot() + geom_line(aes(x = s_start, y = hit_score)) +
  ylim(0, 20) +
  facet_wrap(~CHROM, scales = "free_x")

# nucs <- c("a", "t", "c", "g")
# n_seqs <- 10
# set.seed(313)
# seqs <- matrix(sample(nucs, 20 * n_seqs, replace = T), ncol = 20)
# random_indices <- apply(seqs, MARGIN = 1, function(x) c(1:20)[gRNA_seq_split == x])
# 
# rle(random_indices[[10]])
# random_scores <- lapply(align_indices, function(x) get_score(x))



