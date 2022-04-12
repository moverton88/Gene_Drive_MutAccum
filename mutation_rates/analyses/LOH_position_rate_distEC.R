# We have a set of LOH events and a false homozygous rate, and so we expect
# some LOH events to be false positives, but we do not know which ones.
# In order to perform position-specific analyses without knowing which
# (singleton) LOH events to exlude, we will calculate an LOH/breakpoint
# rate distribution by simulating a random draw of singleton LOH events
# to exclude for each clone given the number of expected false LOH events
all_GT_bounds_merge 
all_LOHbounds_merge

F_LOHrate_ID <- all_GT_bounds_merge %>% 
  group_by(ID) %>% 
  summarize(n_markers = sum(length), 
            n_F_LOH = sum(length) * overall_F_Hom_rate)

n_markers_conv <- all_LOHbounds_merge %>% 
  group_by(ID) %>% 
  summarize(n_markers_conv = sum(length))

marker_stats_xID <- merge(F_LOHrate_ID, n_markers_conv, by = "ID", all = T)
marker_stats_xID$f_markers_conv = marker_stats_xID$n_markers_conv / 
  marker_stats_xID$n_markers

marker_stats_xID %>%
  ungroup() %>% 
  summarize(mean_n_markers = mean(n_markers), 
            min_n_markers = min(n_markers), 
            max_n_markers = max(n_markers),
            mean_n_conv = mean(n_markers_conv), 
            min_n_conv = min(n_markers_conv), 
            max_n_conv = max(n_markers_conv))


loh_cols <- c("Tx_name", "ID", "CHROM", "GT", "start_POSi", 
              "est_start", "est_end", "est_length", 
              "is_error", "isTerm")

singleton_LOH <- all_LOHbounds_merge %>% 
  filter(length == 1) %>% select(any_of(loh_cols))

n_singleton_LOH <- all_LOHbounds_merge %>% 
  filter(length == 1) %>% 
  count(ID)

miss_ID <- evo_IDs[!evo_IDs %in% n_singleton_LOH$ID] 
miss_df <- data.frame(ID = miss_ID, n = 0)
n_singleton_LOH <- rbind(n_singleton_LOH, miss_df)

marker_stats_xID <- merge(marker_stats_xID, n_singleton_LOH, by = "ID")
marker_stats_xID <- marker_stats_xID %>% mutate(n_singles = n) %>% select(!n)

marker_stats_xID %>%
  ggplot() + geom_point(aes(x = n_markers, y = n_markers_conv))

marker_nSingle_lm <- lm(n_singles ~ n_markers, data = marker_stats_xID)
summary(marker_nSingle_lm)

marker_stats_xID %>%
  ggplot() + geom_jitter(aes(x = log10(n_markers), y = n_singles), height = 0.3, width = 0)
  
marker_stats_xID %>%
  ggplot() + geom_jitter(aes(x = n_singles, y = n_markers_conv), height = 0, width = 0.3)


n_singleton_LOH %>% 
  ungroup() %>% 
  summarize(median_sing = median(n), 
            mean_sing = mean(n),
            min_sing = min(n), max_sing = max(n))

marker_stats_xID %>% 
  ungroup() %>% 
  summarize(median_sing = median(n_F_LOH), 
            mean_sing = mean(n_F_LOH),
            min_sing = min(n_F_LOH), max_sing = max(n_F_LOH))



n_singleton_LOH %>% 
  ggplot() + geom_histogram(aes(x = n), binwidth = 1) +
  xlim(-1, 20)

# marker_stats_xID %>% 
#   ggplot() + geom_histogram(aes(x = round(n_F_LOH)), binwidth = 1) +
#   xlim(-1, 20)

data.frame(ID = 1:10000, n = rpois(100, 2.367713)) %>% 
  ggplot() + geom_histogram(aes(x = n), binwidth = 1) +
  xlim(-1, 20)

pnbinom(15, 1, mu = 2.36, lower.tail = F)

data.frame(ID = 1:10000, n = rnbinom(10000, 10, mu = 2.36)) %>% 
  ggplot() + geom_histogram(aes(x = n), binwidth = 1) +
  xlim(-1, 20)

singleton_LOH %>% 
  # group_by(ID) %>%
  summarize(f_BY = sum(GT == "0/0")/n())

multimarker_LOH <- all_LOHbounds_merge %>% 
  filter(length > 1) %>% select(any_of(loh_cols))

n_multimarker_LOH <- all_LOHbounds_merge %>% 
  filter(length > 1) %>% 
  count(ID)

# miss_ID <- evo_IDs[!evo_IDs %in% n_multimarker_LOH$ID] 
# miss_df <- data.frame(ID = miss_ID, n = 0)
# n_multimarker_LOH <- rbind(n_multimarker_LOH, miss_df)

LOHbounds_to_sample <- all_LOHbounds_merge %>% 
  select(ID, length, start_POSi, end_POSi, est_start, est_end)


n_combos <- factorial(n_sLOH)/(factorial(n_err) * factorial(n_sLOH - n_err))
n_perms <- 200

LOHbounds_to_sample$length == 1

# The set of multiple-marker LOHs does not change, so we can split 
# The data into multi and singletons



# One way to measure LOH rate for each clone is to calclulate
# The mean rate at which each marker is converted across simulations

# Another way is to 