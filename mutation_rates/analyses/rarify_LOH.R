

# Rarifaction curve

all_LOHbounds_merge_NS <- all_LOHbounds_merge_NS %>% 
  group_by(ID) %>% 
  mutate(tract_ID = paste0(ID, "_", 1:n()))

tracts_for_markers <- function(SNP_df, LOH_merge_df) {
  # LOH_merge_df <- all_LOHbounds_merge_NS
  tract_ID <- rep("het", nrow(SNP_df))
  for(i in 1:nrow(LOH_merge_df)) {
    s_POSi <- LOH_merge_df$est_start[i]
    e_POSi <- LOH_merge_df$est_end[i]
    id <- as.character(LOH_merge_df$ID[i])
    i_tract <- SNP_df$POSi >= s_POSi & SNP_df$POSi < e_POSi & SNP_df$ID == id
    tract_ID[i_tract] <- LOH_merge_df$tract_ID[i]
  }
  return(tract_ID)
}


tract_IDs <- all_LOHbounds_merge_NS %>% distinct(tract_ID) %>% pull(tract_ID)

LOH_SNPs_ep <- LOH_SNPs %>% filter(Rep != "00")

Line_SNPs <- LOH_SNPs %>% 
  group_by(Line) %>% 
  distinct(POSi, .keep_all = T) %>% 
  arrange(Line, POSi) %>% 
  select(Line, POSi)

LOH_SNPs_ep$tract_ID <- tracts_for_markers(LOH_SNPs_ep, all_LOHbounds_merge_NS)
line_levels <- unique(Line_SNPs$Line)

n_out <- c(100, 500, 1000, 5000, 10000, 15000, 20000, 25000, 30000, 35000, 37000, 38000, 39000, 39500, 40000)
n_reps <- 20
df_out <- data.frame(NULL)
for(l in line_levels) {
  # l <- "F_C"
  LOH_SNPs_ep_line <- LOH_SNPs_ep %>% filter(Line == l)
  line_SNPs <- Line_SNPs %>% filter(Line == l)
  n_markers_line <- nrow(line_SNPs)
  n_markers_LOH <-  LOH_SNPs_ep_line %>% ungroup() %>% filter(tract_ID != "het") %>% nrow()
  line_NULL_df <- data.frame(Line = l,
                             n_anc_m_out = 0, n_anc_m_in = n_markers_line,
                             n_markers_in = nrow(LOH_SNPs_ep_line), sd_n_m_in = 0, 
                             f_markers_in = 1, 
                             n_markers_LOH_in = n_markers_LOH, f_markers_LOH_in = n_markers_LOH/nrow(LOH_SNPs_ep_line),
                             n_LOH_in =  length(unique(LOH_SNPs_ep_line$tract_ID)),
                             sd_n_LOH_in = 0)
  line_df <- data.frame(NULL)
  for(j in 1:length(n_out)){
    # j <- 2
    if(n_out[j] < n_markers_line) {
      # sample markers to remove across end-point clones
      posi_out_matrix <- replicate(n_reps, sample(line_SNPs$POSi, n_out[j]))
      
      n_markers_in <- c()
      n_markers_LOH_in <- c()
      n_LOH_in <- c()
      for(i in 1:n_reps) {
        # i = 1
        n_markers_in[i] <- LOH_SNPs_ep_line %>% 
          filter(POSi %in% posi_out_matrix[, i]) %>% 
          nrow()
        
        n_markers_LOH <- LOH_SNPs_ep_line %>% 
          filter(POSi %in% posi_out_matrix[, i], 
                 substr(tract_ID, nchar(tract_ID) - 2, nchar(tract_ID)) != "het") %>% 
          count(tract_ID) %>% filter(n > 1)
        
        n_markers_LOH_in[i] <- sum(n_markers_LOH$n)
        
        n_LOH_in[i] <- n_markers_LOH %>%
          distinct(tract_ID) %>% 
          nrow()
        
        if(i == n_reps) {
          mean_n_df <- data.frame(Line = l, 
                                  n_anc_m_out = n_out[j], n_anc_m_in = n_markers_line - n_out[j],
                                  n_markers_in = mean(n_markers_in), sd_n_m_in = sd(n_markers_in), 
                                  f_markers_in = mean(n_markers_in)/line_NULL_df$n_markers_in, 
                                  n_markers_LOH_in = mean(n_markers_LOH_in), 
                                  f_markers_LOH_in = mean(n_markers_LOH_in/n_markers_in),
                                  n_LOH_in = mean(n_LOH_in),
                                  sd_n_LOH_in = sd(n_LOH_in))
        }
      }
      line_df <- rbind(line_df, mean_n_df)
    }
  }
  line_df <- rbind(line_NULL_df, line_df)
  df_out <- rbind(df_out, line_df)
  print(l)
}

rarified_lines <- df_out

rarified_all <- df_out %>% group_by(n_anc_m_out) %>% summarize(tot_LOH_in = sum(n_LOH_in))

# line_df_1 <- line_df
df_out_line <- df_out %>% filter(Line == "N_C") 

X <- rarified_all %>%
  # filter(n_anc_m_out < 30000, n_anc_m_out > 0) %>% 
  pull(n_anc_m_out)

a <- rarified_all$tot_LOH_in[1]
a <- 1200
b <- 0
c <- 0.0003
d <- 8.3E-5*9/560

rarified_all$pred <- (a - (a - b) * exp (-c * X)) + d * X

rarified_all %>% 
  mutate(n_LOH_in = tot_LOH_in) %>%
  filter(n_anc_m_out < 30000, n_anc_m_out > 0) %>%
  # group_by(Line) %>%
  # mutate(f_LOH_in = n_LOH_in/max(n_LOH_in)) %>%
  # filter(f_markers_in != 1) %>%
  ggplot(aes(x = log10(n_anc_m_out), y = n_LOH_in)) + 
  # ggplot(aes(x = n_anc_m_out, y = n_LOH_in, color = Line)) + 
  geom_point() + 
  # geom_errorbar(aes(ymin = n_LOH_in - sd_n_LOH_in, ymax = n_LOH_in + sd_n_LOH_in)) +
  geom_line() #+
  # geom_line(aes(y = pred), color = 'blue')

rarified_lines %>% 
  # mutate(n_LOH_in = tot_LOH_in) %>%
  filter(n_anc_m_in > 10^3.5, n_anc_m_out > 0) %>%
  # group_by(Line) %>%
  # mutate(f_LOH_in = n_LOH_in/max(n_LOH_in)) %>%
  # filter(f_markers_in != 1) %>%
  # ggplot(aes(x = log10(n_anc_m_in), y = n_LOH_in)) + 
  ggplot(aes(x = log10(n_anc_m_in), y = n_LOH_in, color = Line)) +
  geom_point() + 
  # geom_errorbar(aes(ymin = n_LOH_in - sd_n_LOH_in, ymax = n_LOH_in + sd_n_LOH_in)) +
  geom_line() #+
geom_line(aes(y = pred), color = 'blue')

(4, 20) 4.45, 30
10/0.45
