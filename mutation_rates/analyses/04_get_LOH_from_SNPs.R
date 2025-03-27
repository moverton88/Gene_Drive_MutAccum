


final_cols <- c("CHROM", "POS", "POSi", "Tx_name", "Tx_ID", "Tx", "Line", "Rep", "ID", 
                "GT", "GQ", "Ref_DP_final", "Alt_DP_final", "Sum_DP_final",
                "existing_SNP", "Cut")

# For each founder group, collect calls in all clones that are Het in the founder
LOH_SNPs <- SNPs_merge_finalGT %>% 
  ungroup() %>%
  filter(!POSi %in% high_error_POSi, 
         QUAL_BYcall >= 1000, QUAL_RMcall >= 1000, 
         !Cut) %>% 
  select(all_of(final_cols)) %>% 
  anc_GT_fltr(anc_GT = "0/1")
LOH_SNPs$Line <- droplevels(LOH_SNPs$Line)
LOH_SNPs$ID <- droplevels(LOH_SNPs$ID)



# Tables of LOH boundaries and counts for each clone ------
all_GT_bounds <- LOH_SNPs %>% 
  filter(Rep != "00") %>% 
  filter(POSi %in% clean_markers) %>%
  EstDataBounds(., rm_noData = T)

all_GT_bounds_aneuCorr <- all_GT_bounds %>% 
  mutate(ID_CHROM = paste0(ID, "_", CHROM)) %>% 
  filter(!ID_CHROM %in% ID_CHROM_aneuploid$ID_CHROM)

# Mark LOH regions belonging to complex events
all_GT_bounds_aneuCorr <- all_GT_bounds_aneuCorr %>% 
  MarkLOHcomplexes(., gap = 10000)

all_GT_bounds_merge <- all_GT_bounds_aneuCorr %>%
  MergeComplexLOHs() %>% 
  MarkTerminalLOHs(ancHet = LOH_SNPs)

save(all_LOHbounds_merge_NS, file = bounds_filename)



