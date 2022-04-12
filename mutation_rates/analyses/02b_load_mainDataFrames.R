# Loads in tables generated in 02_mainDataFrames.R


# Table of merged callsets with final genotyping (SNPs_merge_finalGT)
###############################################################################

SNPs_ccs <- c("NULL", "factor", "integer", "integer", "factor", 
         "factor", "factor", "factor", 
         "numeric", "character", "character", "factor", "integer", "integer", "integer", 
         "numeric", "character", "character", "factor", "integer", "integer", "integer", 
         "logical", "factor", "integer", "integer", "integer", "integer", 
         "factor", "logical", "numeric", "logical")

SNPs_merge_finalGT <- read.delim(SNPs_merge_filename, header = T, sep = " ", quote = "\"",
                         colClasses = SNPs_ccs, row.names = NULL)


SNPs_merge_finalGT$Tx <- factor(SNPs_merge_finalGT$Tx, levels = Tx_levels)
SNPs_merge_finalGT$Tx_name <- factor(SNPs_merge_finalGT$Tx_name, levels = Tx_name_levels)

make_ref_calls_final <- function(df, ref = "BY") {
  # df <- head(SNPs_merge_finalGT)
  col_IDs <- c("GT", "GQ", "Ref_DP", "Alt_DP")
  ref_cols <- paste0(col_IDs, "_", ref, "call")
  final_cols <- c("GT", "GQ", "Ref_DP_final", "Alt_DP_final")
  df[, final_cols] <- df[, ref_cols]
  df$Sum_DP_final <- df$Ref_DP_final + df$Alt_DP_final
  return(df)
}

# SNPs_merge_finalGT <- SNPs_merge_finalGT %>% make_ref_calls_final()


all_error_rates <- SNPs_merge_finalGT %>% 
  filter(GQ >= 50, !Cut, !Line %in% noAncestor, QUAL_BYcall >= 1000, QUAL_RMcall >= 1000) %>%
  errorFromPhylo(flsHom_support = 4, flsHet_support = 6, output_POSi = F)

overall_F_Hom_rate <- sum(all_error_rates$n_F_Hom)/(sum(all_error_rates$n_Het_q))
overall_F_Het_rate <- sum(all_error_rates$n_F_Het)/sum(all_error_rates$n_Hom_q)
overall_F_Hom_rate
overall_F_Het_rate

all_LOHbounds <- LOH_SNPs %>% 
  filter(Rep != "00") %>% 
  filter(POSi %in% clean_markers) %>%
  EstDataBounds(., rm_noData = T)

# Mark LOH regions belonging to complex events
all_LOHbounds <- all_LOHbounds %>% 
  MarkLOHcomplexes(., gap = 10000) 

all_LOHbounds_merge <- all_LOHbounds %>%
  MergeComplexLOHs() 

all_LOHbounds_merge_EC <- all_LOHbounds_merge %>% 
  MarkErrorLOH(., error_rate = overall_F_Hom_rate)

# all_LOHbounds_merge_EC <- all_LOHbounds_merge_EC %>% 
#   select(-is_error) %>%
#   MarkErrorLOH(., error_rate = overall_F_Hom_rate)

all_LOHbounds_merge_EC <- all_LOHbounds_merge_EC %>% 
  MarkTerminalLOHs(ancHet = LOH_SNPs)

all_LOHcounts_merge_EC <- all_LOHbounds_merge_EC %>% 
  CountLOHevents(omitError = T)
all_LOHcounts_merge_EC <- CategoriesFromID(all_LOHcounts_merge_EC_F_D)

all_LOHcounts_merge <- all_LOHbounds_merge_EC %>% 
  CountLOHevents(omitError = F)



SNPs_merge_finalGT %>% filter(ID == "H_C03", POSi > 5392574, POSi < 5824339, GT == "0/1")
all_LOHbounds_merge_EC %>% filter(ID == "H_C03", GT != "0/1")
H_C_anc_SNPs <- SNPs_merge_finalGT %>% filter(ID == "H_C00", GT == "0/1", !Cut) %>% pull(POSi)

SNPs_merge_finalGT %>% 
  filter(ID == "H_D03") %>%
  # filter(Tx == "H") %>%
  filter(POSi %in% H_D_anc_SNPs) %>%
  group_by(CHROM) %>% summarize(f_alt = median(f_Alt)) %>%
  # filter(CHROM == "07") %>%
  # filter(POSi > 800000, POSi < 900000) %>%
  # filter(GT == "0/0") %>%
  ggplot(aes(x = CHROM, y = f_alt)) + geom_line() #+ 
  facet_grid(ID~.)

SNPs_merge_finalGT %>% filter(ID == "H_C03", GT == "0/0") %>% ggplot() + geom_histogram(aes(x = GQ))
SNPs_merge_finalGT %>% 
  # filter(ID == "H_C03") %>% 
  group_by(GT) %>% summarize(m = mean(GQ))

sitewise_GTs <- SNPs_merge_finalGT %>% 
  filter(GQ >= 50, !repeats, !Cut) %>%
  site_genotype_stats(group = "all")
  
LOH_poor_GT_sites <- sitewise_GTs %>% filter(fHet_anc < 0.9) %>% pull(POSi)

# Table of SNPs used for LOH analysis
###############################################################################

LOH_ccs <- c("NULL", "factor", "integer", "integer", "factor", 
             "factor", "factor", "factor", "factor", "factor", 
             "numeric", "logical", "logical", "integer")

LOH_SNPs <- read.delim(LOH_SNPs_file, header = T, sep = " ", quote = "\"",
                                 colClasses = LOH_ccs, row.names = NULL)

LOH_SNPs$Tx <- factor(LOH_SNPs$Tx, levels = Tx_levels)
LOH_SNPs$Tx_name <- factor(LOH_SNPs$Tx_name, levels = Tx_name_levels)
n_clones_xTx <- LOH_SNPs %>% distinct(ID, .keep_all = T) %>% count(Tx_name)
n_evo_xTx <- LOH_SNPs %>% filter(Rep != "00") %>% distinct(ID, .keep_all = T) %>% count(Tx_name)
# Table of de novo SNM mutations
###############################################################################

DN_ccs <- c("NULL", "factor", "integer", "integer", "factor", 
             "factor", "factor", "factor", "factor", "factor", 
             "numeric", "logical", "logical")

denovo_SNPs <- read.delim(denovo_SNPs_file, header = T, sep = " ", quote = "\"",
                       colClasses = DN_ccs, row.names = NULL)

denovo_SNPs$Tx <- factor(denovo_SNPs$Tx, levels = Tx_levels)
denovo_SNPs$Tx_name <- factor(denovo_SNPs$Tx_name, levels = Tx_name_levels)

# Table of position bounds of regions with contiguous genotypes for LOH analysis
###############################################################################

bounds_ccs <- c("NULL", "integer", "integer", "factor", "integer", 
                "integer", "integer", "integer", "integer", "integer", 
                "integer", "factor", "factor", "factor", "factor", 
                "factor", "factor", "integer", "integer", "logical", 
                "logical")

all_LOHbounds_merge_EC <- read.delim(bounds_filename, header = T, sep = " ", quote = "\"",
                          colClasses = bounds_ccs, row.names = NULL)

all_LOHbounds_merge_EC$Tx <- factor(all_LOHbounds_merge_EC$Tx, levels = Tx_levels)
all_LOHbounds_merge_EC$Tx_name <- factor(all_LOHbounds_merge_EC$Tx_name, levels = Tx_name_levels)

# Table of position bounds of regions with contiguous genotypes for LOH analysis
###############################################################################
counts_ccs <- c("NULL", "factor", "integer", "integer", "integer", 
                "integer", "integer", "integer", "integer", "integer", 
                "integer", "integer", "factor", "factor", "factor", 
                "factor")

all_LOHcounts_merge_EC <- read.delim(countsEC_filename, header = T, sep = " ", quote = "\"",
                                     colClasses = counts_ccs, row.names = NULL)

all_LOHcounts_merge_EC$Tx <- factor(all_LOHcounts_merge_EC$Tx, levels = Tx_levels)
all_LOHcounts_merge_EC$Tx_name <- factor(all_LOHcounts_merge_EC$Tx_name, levels = Tx_name_levels)

all_LOHcounts_merge <- read.delim(counts_filename, header = T, sep = " ", quote = "\"",
                               colClasses = counts_ccs, row.names = NULL)

all_LOHcounts_merge$Tx <- factor(all_LOHcounts_merge$Tx, levels = Tx_levels)
all_LOHcounts_merge$Tx_name <- factor(all_LOHcounts_merge$Tx_name, levels = Tx_name_levels)


