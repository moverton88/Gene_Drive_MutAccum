
###############################################################################

ColSplit <- function(col, pat, nms) {
  # col <- c("youbad", "youRgood", "youRaverage", "youRnotRgoingRtoRwork")
  # pat <- "R"
  # nms <- c("first", "second")
  col_list <- sapply(col, function(x) strsplit(x, pat))
  col_list_fill <- sapply(col_list, function(x) 
                                      ifelse(length(x) >= length(nms), 
                                             list(x[1:length(nms)]), 
                                             list(c(x, rep(NA, length(nms) - length(x))))))
  cols_out <- data.frame(t(data.frame(col_list_fill)), row.names = NULL, stringsAsFactors = F)
  colnames(cols_out) <- nms
  return(cols_out)
}

OperateByFactorMatch <- function(op_pair, trgt_pair, .fun) {
  # Each input is a pair of columns to operate over, with the first of 
  # each pair the levels and the second the operator and operand. 
  # .fun is a function of the form function(x, y) {x + y}
  
  # op_pair <- chrom_bounds_Ref %>% select(CHROM, Start)
  # trgt_pair <- all_var_P1 %>% mutate(diff_POSs = POSi - POS) %>% select(CHROM, diff_POSs)
  # op_pair <- as.data.frame(op_pair)
  # trgt_pair <- as.data.frame(trgt_pair)
  op_levels <- levels(op_pair[, 1])
  trgt_levels <- levels(trgt_pair[, 1])
  trgt_pair$rn <- 1:nrow(trgt_pair)
  trgt_names <- colnames(trgt_pair)
  if(length(op_levels) == length(trgt_levels)) {
    if(all(op_levels %in% trgt_levels)) {
      op_pair[, 1] <- factor(op_pair[, 1], trgt_levels)
      trgt_out <- data.frame(NULL)
      for(lvl in trgt_levels) {
        # lvl <- trgt_levels[1]
        trgt_sub_in <- trgt_pair[trgt_pair[1] == lvl, 2:3]
        op_sub <- op_pair[op_pair[1] == lvl, 2]
        # .fun <- function(x, y) {x + y}
        trgt_sub_out <- .fun(op_sub, trgt_sub_in[, 1])
        trgt_sub_out <- data.frame(value = trgt_sub_out, trgt_sub_in[2])
        trgt_out <- rbind(trgt_out, trgt_sub_out)
      }
      trgt_out <- trgt_out[order(trgt_out[, 2]), 1]
      colnames(trgt_pair) <- trgt_names
      return(trgt_out)
    } else {
      print("factor levels do not match")
    }
  } else {
    print("factor levels are different lengths")
  }
}

GetChromDF <- function(chrom_lengths) {
  # chrom_lengths <- chr_info %>% 
  #   select(Qry_CHROM, Qry_length) %>% 
  #   rename(CHROM_acc = Qry_CHROM, length = Qry_length) %>%
  #   mutate(CHROM = str_pad(1:nrow(chr_info), width = 2, pad = "0"))
  chrom_lengths$CHROM <- factor(chrom_lengths$CHROM)
  n_chroms <- length(levels(chrom_lengths$CHROM))
  chrom_num <- as.numeric(chrom_lengths$CHROM)
  roman_chr <- factor(as.character(as.roman(chrom_num)), 
                      levels = as.character(as.roman(chrom_num)))
  chrom_IDs <- data.frame(CHROM = chrom_lengths$CHROM, 
                          rom_CHROM = roman_chr)
  chrom_lengths <- merge(chrom_IDs, chrom_lengths, by = "CHROM")
  
  # chrom_indcs <- cumsum(chrom_lengths$length)
  
  chr_start <- cumsum(c(1, chrom_lengths$length[1:(n_chroms - 1)]))
  chr_end <- cumsum(chrom_lengths$length)
  
  chr_bound <- data.frame(CHROM = chrom_lengths$CHROM, 
                          Start = as.numeric(chr_start), 
                          End = as.numeric(chr_end))
  
  # chrom_bounds <- BoundsFromLengths(chrom_lengths$length)
  # chrom_bounds <- cbind(chrom_lengths, chrom_bounds[c("Start", "End")])
  chrom_bounds <- merge(chrom_lengths, chr_bound, by = "CHROM")
  return(chrom_bounds)
}

ChainToDFlist <- function(chain_file, n_chrs = 16)  {
  ###########################################################################
  # A chain file is composed of a series of sections, one for each contig or 
  # chromosome. Within each section, there are three columns:
  # 1) block length, 
  # 2) difference between end of current block and start of next block in the 
  # reference sequence, 
  # 3) difference between end of this bloc and beginning of next query block
  # This function reads in the chain file and converts it to a data.frame 
  # capable of converting the position coordinates of the query vcf back to 
  # that of the Type reference sequence.
  # It also creates data.frames of the chromosome lengths for both the
  # reference and query strains from the chain file.
  # The function outputs a list of these three data.frames
  
  # chain_file <- P1_chainFile
  chain_table <- readLines(chain_file) # read in chain file
  chain_table <- gsub(pattern = "\"", replacement = "", x = chain_table) # remove spurious "\" characters
  chain_table <- chain_table[!chain_table == ""] # remove new lines at end of each section
  chr_idx <- grep("chain", chain_table) # index section headers
  if(length(chr_idx) > n_chrs) {
    # remove extra chromosomes if present
    last_idx <- chr_idx[(n_chrs + 1)]
    chain_table <- chain_table[1:(last_idx - 1)]
    chr_idx <- chr_idx[1:n_chrs]
  } 
  end_idx <- length(chain_table)
  value_idx <- data.frame(from = chr_idx + 1, 
                          to = c(chr_idx[2:length(chr_idx)] - 1, end_idx)) # dataframe of section blocks
  chr_headers <- chain_table[chr_idx] # get section headers
  chr_info <- str_split(chr_headers, " ", simplify = T) %>% as.data.frame(stringsAsFactors = F) # turn headers into dataframe
  col_names <- c("CHROM", "length", "strand", "start", "end") # rename dataframe colnames
  colnames(chr_info) <- c("null", "score", paste0("Ref_", col_names), paste0("Qry_", col_names), "ID")
  num_cols <- c(paste0("Ref_", col_names[c(2, 4:5)]), paste0("Qry_", col_names[c(2, 4:5)]))
  chr_info[, num_cols] <- sapply(chr_info[, num_cols], as.numeric)
  chrom_names <- str_pad(1:nrow(chr_info), width = 2, side = "left", pad = "0") # get chromosome names
  
  Ref_chrom_starts <- cumsum(c(1, chr_info$Ref_length[1:(nrow(chr_info) - 1)]))
  Qry_chrom_starts <- cumsum(c(1, chr_info$Qry_length[1:(nrow(chr_info) - 1)]))
  
  chain_df <- data.frame(NULL)
  for (ch in 1:nrow(chr_info)) {
    # ch = 1
    # subName <- chain_table[chr_idx[ch]]
    sub_chain <- chain_table[value_idx$from[ch]:value_idx$to[ch]]
    sub_df <- data.frame(sapply(ColSplit(sub_chain, " ", c("block", "Ref_diff", "Qry_diff")), as.numeric))
    sub_df[is.na(sub_df)] <- 0
    sub_df$cs_block <- cumsum(sub_df$block)
    sub_df$cs_Ref_diff <- cumsum(sub_df$Qry_diff) - cumsum(sub_df$Ref_diff)
    sub_df$cs_Qry_diff <- cumsum(sub_df$Ref_diff) - cumsum(sub_df$Qry_diff)
    sub_df$Ref_POS <- sub_df$cs_block + cumsum(sub_df$Ref_diff) + 1 # cs_block is 0 indexed, so add 1 for position
    sub_df$Qry_POS <- sub_df$cs_block + cumsum(sub_df$Qry_diff) + 1
    # sub_df$Ref_POS[1] <- sub_df$cs_block[1] + 1
    # sub_df$Qry_POS[1] <- sub_df$cs_block[1] + 1
    sub_df[1, c("Ref_POS", "Qry_POS")] <- sub_df$cs_block[1] + 1
    sub_df$Ref_POSi <- sub_df$Ref_POS + Ref_chrom_starts[ch] - 1 # chrom_starts is 1 indexed, so -1 to add intervening chroms
    sub_df$Qry_POSi <- sub_df$Qry_POS + Qry_chrom_starts[ch] - 1
    sub_df$CHROM <- chrom_names[ch]
    chain_df <- rbind(chain_df, sub_df)
  }
  chain_df$CHROM <- factor(chain_df$CHROM)
  
  chrom_info_Ref <- chr_info %>% 
    select(Ref_CHROM, Ref_length) %>% 
    rename(CHROM_acc = Ref_CHROM, length = Ref_length) %>%
    mutate(CHROM = str_pad(1:nrow(chr_info), width = 2, pad = "0")) %>% 
    GetChromDF()
  
  chrom_info_Qry <- chr_info %>% 
    select(Qry_CHROM, Qry_length) %>% 
    rename(CHROM_acc = Qry_CHROM, length = Qry_length) %>%
    mutate(CHROM = str_pad(1:nrow(chr_info), width = 2, pad = "0")) %>% 
    GetChromDF()
  
  chain_list <- list(chain = chain_df, Ref_chrom = chrom_info_Ref, Qry_chrom = chrom_info_Qry)
  
  return(chain_list)
}

ConvertPosIndicies <- function(pos_df, pos_col = "POS", chrom_col = "CHROM", 
                               index_out = "POSi", chrom_bound,
                               add_chroms = F) {
  # pos_df <- POSi_data_in
  # pos_col = "est_start"
  # index_out = "POS"
  pos_df <- as.data.frame(pos_df)
  pos_out <- pos_df[, pos_col]
  chrom_out <- rep("00", nrow(pos_df))
  if(index_out == "POSi") {
    for(ch in 1:nrow(chrom_bound)) {
      chrom <- chrom_bound$CHROM[ch]
      s <- chrom_bound$Start[ch]
      i_sub <- pos_df$CHROM == chrom
      pos_out[i_sub] <- pos_df[i_sub, pos_col] + s - 1
    }
  }
  if(index_out == "POS") {
    for(ch in 1:nrow(chrom_bound)) {
      # ch = 5
      s <- chrom_bound$Start[ch]
      e <- chrom_bound$End[ch]
      i_sub <- pos_df[, pos_col] >= s & pos_df[, pos_col] <= e
      pos_out[i_sub] <- pos_df[i_sub, pos_col] - s + 1
      if(add_chroms) {
        chrom_out[i_sub] <- as.character(chrom_bound$CHROM[ch])
      }
    }
  }
  if(add_chroms) {
    list_out <- list(pos_out, chrom_out)
    names(list_out) <- c(index_out, "CHROM")
    return(list_out)
  } else {
    return(pos_out)
  }
}

HC_multiVcfParse_allVar <- function(vcf_file., max_chrs = 16, verbose = F) {
  ## takes as input a VCF file output from GATK CombineGVCFs => GenotypeGVCFs
  # This file contains all samples from a founder group, with each sample's genotype
  # information occupying a single column
  # vcf_file. <- P1_file
  vcf_in <- read.table(vcf_file., stringsAsFactors = F)
  cat("VCF read in\n")

  # Extract column header from commented section of VCF
  hdr_full <- readLines(vcf_file., n = 120)
  skip <- max(grep("^\\s*##", hdr_full))
  # Extract header line, tab-separated values to vector,
  # remove "#" from CHROM, name as dataframe columns
  hdr_raw <- hdr_full[(skip + 1)]
  # Split header text into column name vector
  hdr_vctr <- unlist(strsplit(hdr_raw, split="\t"))
  # Remove "#" from CHROM
  hdr_vctr[1] <- substr(hdr_vctr[[1]], start = 2, stop = 6)
  # Set column names
  colnames(vcf_in) <- hdr_vctr
  # Remove unneeded columns
  vcf_in <- select(vcf_in, -ID, -FILTER, -INFO, -FORMAT)
  
  # Get chromosome length info
  contigs <- hdr_full[grep("^\\s*##contig=", hdr_full)]
  contig_lengths = sapply(strsplit(contigs, ","), "[[", 2)
  chrom_lengths <- as.numeric(gsub("[^0-9]", "", contig_lengths))[1:max_chrs]
  vcf_in$chrom_n <- as.numeric(factor(vcf_in$CHROM))
  vcf_in$CHROM <- factor(str_pad(vcf_in$chrom_n, 2, pad = "0"))

  if(max(vcf_in$chrom_n) >= (max_chrs + 1)) {
    vcf_in <- filter(vcf_in, !CHROM %in% c(as.character((max_chrs + 1):max(vcf_in$chrom_n))))
  }

  chrom_indcs <- cumsum(c(0, chrom_lengths[1:(max_chrs - 1)]))
  vcf_in$POSi <- vcf_in$POS + chrom_indcs[as.numeric(vcf_in$CHROM)]

  sampleColNames <- c("GT", "AlleleDP", "Sum_DP", "GQ", "PGT", "PID", "PhredLike", "PS")
  nonsample_cols <- c("CHROM", "POS", "REF", "ALT", "QUAL", "FORMAT", "POSi")
  vcf_long <- vcf_in %>% 
                pivot_longer(cols = !contains(nonsample_cols), 
                             names_to = "ID", values_to = "raw_sample_info") %>%
                arrange(ID, POSi)

  vcf_long <- vcf_long %>%
              mutate(raw_sample_info = gsub("|", "/", raw_sample_info, fixed = TRUE)) %>%
              separate(raw_sample_info, into = sampleColNames, sep = ":", 
                        extra = "merge", fill = "right") %>%
              separate(ALT, into = c("ALT", "ALT_excess"), sep = ",", 
                        extra = "merge", fill = "right") %>%
              separate(AlleleDP, into = c("Ref_DP", "Alt_DP", "DP_excess"), sep = ",", 
                        extra = "merge", fill = "right") %>%
              mutate(across(c(Ref_DP, Alt_DP), as.numeric),
                      Sum_DP = Ref_DP + Alt_DP) %>%
              select(-contains("excess"), -contains(sampleColNames[5:8]))

  vcf_long <- cbind(select(vcf_long, CHROM:POS),
                  select(vcf_long, POSi),
                  select(vcf_long, REF:QUAL),
                  select(vcf_long, GT:GQ),
                  select(vcf_long, ID),
                  stringsAsFactors = FALSE)

  vcf_long <- filter(vcf_long, GT %in% c("0/0", "0/1", "1/1"))
  vcf_long$GT <- factor(vcf_long$GT, levels = c("0/0", "0/1", "1/1"))
  vcf_long$GQ <- as.numeric(vcf_long$GQ)
  vcf_long$CHROM <- factor(vcf_long$CHROM)
  return(vcf_long)
}


# There should be another version that gets the chrom lengths from the chain file
ConstructLiftover <- function(vcf_df, chain_df,
                           chrom_from, chrom_to,
                           pos_col = "Qry_POS",
                           csDiff_col = "cs_Qry_diff") {
  # The chain_df, chrom_from, and chrom_to, arguements 
  # should be the three elements in the ChainToDFlist output in order
  # vcf_df <- all_var_P1
  # chain_df <- P1toRef_chainDf
  # pos_col <- "Qry_POS"
  # csDiff_col <- "cs_Qry_diff"
  # chrom_from = chrom_bounds_P1
  # chrom_to = chrom_bounds_Ref
  vcf_df$POS_adj <- vcf_df$POS
  vcf_df$POSi_adj <- vcf_df$POSi
  # chain_df$gDiff <- chain_df[, gDiff_col_1] - chain_df[, gDiff_col_2]
  for(ch in levels(chain_df$CHROM)) {
    # ch = levels(chain_df$CHROM)[1]
    ch_chain <- chain_df %>% filter(CHROM == ch)
    chDiff <- chrom_to$Start[chrom_to$CHROM == ch] - chrom_from$Start[chrom_from$CHROM == ch]
    # iChrom_vcf <- vcf_df$CHROM == ch
    nBlocks <- nrow(ch_chain)
    for(b in 1:(nBlocks - 1)) {
      # b = 3
      # We isolate positions located within each block
      # Each block is defined as an uninterrupted alignment
      # All positions up to Qry_POSi are aligned, and so the positions
      # between Qry_POSi$cs_Qry_diff[b] and Qry_POSi$cs_Qry_diff[b+1] 
      # need to be adjusted by the value of cs_Qry_diff[b]
      # min_POSi <- ifelse(b == 1, 0, ch_chain[(b - 1), posi_col])
      b_i <- vcf_df$CHROM == ch & 
              vcf_df$POS >= ch_chain[b, pos_col] & 
              vcf_df$POS < ch_chain[(b + 1), pos_col]
      vcf_df$POS_adj[b_i] <- vcf_df$POS[b_i] + ch_chain[b, csDiff_col]
      # vcf_df$POSi_adj[b_i] <- vcf_df$POSi[b_i] + ch_chain[b, csDiff_col] + chDiff
    }
    # e_i <- vcf_df$POSi >= ch_chain[nBlocks, posi_col] & vcf_df$POSi <= chrom_from$End[as.numeric(ch)]
    # vcf_df$POS_adj[e_i] <- vcf_df$POS[e_i] + ch_chain[nBlocks, csDiff_col]
    # lDiff <-  chrom_to$Start[nrow(chrom_to)] - chrom_from$Start[nrow(chrom_to)]
    # vcf_df$POSi_adj[e_i] <- vcf_df$POSi[e_i] + ch_chain[nBlocks, csDiff_col] + lDiff
    print(ch)
  }
  # vcf_df %>% 
  #   mutate(diff_POSs = POSi - POS) %>% 
  #   ggplot(aes(x = POS, y = diff_POSs, color = ID)) + 
  #   geom_point() + facet_wrap(~CHROM, scales = "free")
  
  vcf_df$POSi_adj <- ConvertPosIndicies(vcf_df, pos_col = "POS_adj", chrom_col = "CHROM", 
                                 index_out = "POSi", chrom_bound = chrom_to,
                                 add_chroms = F)
  # Separate data by chromosome and count number of sites for which POS and POSi do not agree
  n_POSs_diff <- OperateByFactorMatch(chrom_to %>% select(CHROM, Start), 
                                  vcf_df %>% mutate(diff_POSs = POSi_adj - POS_adj) %>% select(CHROM, diff_POSs),
                                  .fun = function(x, y) {x - y - 1}) %>% sum(. != 0)
  
  print(paste0(n_POSs_diff, "/", nrow(vcf_df), "sites where POS and POSi do not agree"))
  vcf_df$POS <- vcf_df$POS_adj
  vcf_df$POSi <- vcf_df$POSi_adj
  vcf_df <- vcf_df %>% select(-c(POS_adj, POSi_adj))
  return(vcf_df)
}

SwapCrossedAlleles <- function(SNP_df, annotated_POSi = NULL, to_match = "P1call",
                               to_swap = "P2call") {
  # We expect crossed alleles at sites that differ between the two parental sequences.
  # This may be the same set of sites as is present in a vcf of variants between the 
  # two reference sequences
  # We can also simply swap alleles for parental variants present in the empirical data
  
  # SNP_df <- SNPs_merge
  # annotated_POSi <- RMxBY_comp_SNPs_POSi
  
  # Index columns to match and to swap by searching for keywords
  cnm <- colnames(SNP_df)
  
  # Columns referring to which parental reference
  i_match <- grep(to_match, cnm)
  i_swap <- grep(to_swap, cnm)
  
  # Columns referring to alleles
  i_ref_alel <- grep("REF", cnm)
  i_alt_alel <- grep("ALT", cnm)
  
  # Construct column index for each allele and parent
  i_alel_ord <- c(intersect(i_match, i_ref_alel), 
                  intersect(i_match, i_alt_alel),
                  intersect(i_swap, i_ref_alel),
                  intersect(i_swap, i_alt_alel))
  
  # Columns with allele depth
  i_ref_DP <- grep("Ref_DP", cnm)
  i_alt_DP <- grep("Alt_DP", cnm)
  i_DP_ord <- c(intersect(i_match, i_ref_DP), 
                intersect(i_match, i_alt_DP),
                intersect(i_swap, i_ref_DP),
                intersect(i_swap, i_alt_DP))
  
  # Columns with genotypes
  i_GT <- grep("GT", cnm)
  i_GT_ord <- c(intersect(i_match, i_GT), 
                intersect(i_swap, i_GT))
  
  if(length(c(i_alel_ord, i_DP_ord, i_GT_ord)) != 10) {
    print("Table does not contain the correct columns")
  }
  
  # Index sites where the REF allele of one callset matches the ALT allele of the other
  # and vis-versa
  iCross_alleles <- SNP_df[, i_alel_ord[1]] == SNP_df[, i_alel_ord[4]] & 
    SNP_df[, i_alel_ord[2]] == SNP_df[, i_alel_ord[3]]
  
  # Do not include sites with NAs by assigning as FALSE
  i_na_alleles <- is.na(iCross_alleles)
  iCross_alleles[is.na(iCross_alleles)] <- F
  
  if(length(annotated_POSi) != 0) {
    i_Pvariants <- SNP_df$POSi %in% annotated_POSi
    # i_Pvariants[3:4] <- T
    
    # Record POSi for which the alleles are crossed but missing in the Pvariant df or
    # for which they are present in the Pvariant df but not crossed as expected
    i_dscd_source <- (!i_Pvariants & iCross_alleles) | (i_Pvariants & !iCross_alleles)
  }
  
  # The set of sites to swap alleles: Known from references or data and not 
  # missing in the P2 callset
  iSwap_alleles <- iCross_alleles & !is.na(SNP_df$REF_P2call)
  
  # Make copies of the allele columns to be edited and then perform swap
  SNP_df$REF_swap <- SNP_df[, i_alel_ord[3]]
  SNP_df$ALT_swap <- SNP_df[, i_alel_ord[4]]
  SNP_df[iSwap_alleles, i_alel_ord[4]] <- SNP_df$REF_swap[iSwap_alleles]
  SNP_df[iSwap_alleles, i_alel_ord[3]] <- SNP_df$ALT_swap[iSwap_alleles]
  
  SNP_df$Ref_DP_swap <- SNP_df[, i_DP_ord[3]]
  SNP_df$Alt_DP_swap <- SNP_df[, i_DP_ord[4]]
  SNP_df[iSwap_alleles, i_DP_ord[4]] <- SNP_df$Ref_DP_swap[iSwap_alleles]
  SNP_df[iSwap_alleles, i_DP_ord[3]] <- SNP_df$Alt_DP_swap[iSwap_alleles]
  
  # Copy GT column
  SNP_df$GT_swap <- SNP_df[, i_GT_ord[2]]
  
  # Index homozygous sites in P2 callset
  iAlt_P2call <- SNP_df[, i_GT_ord[2]] == "0/0"
  iRef_P2call <- SNP_df[, i_GT_ord[2]] == "1/1"
  
  # Edit original GT column to swap genotypes of crossed alleles
  SNP_df[iSwap_alleles & iAlt_P2call, i_GT_ord[2]] <- "1/1"
  SNP_df[iSwap_alleles & iRef_P2call, i_GT_ord[2]] <- "0/0"
  
  # Index matching alleles after swap
  iSame_het <- SNP_df[, i_alel_ord[1]] == SNP_df[, i_alel_ord[3]] & 
                SNP_df[, i_alel_ord[2]] == SNP_df[, i_alel_ord[4]]
  
  # Column of sites which still do not match correctly
  
  if(length(annotated_POSi) != 0) {
    SNP_df$bad_alleles <- !iSame_het | i_dscd_source | i_na_alleles
  } else {
    SNP_df$bad_alleles <- !iSame_het | i_na_alleles
  }
  
  SNP_df <- SNP_df %>% select(!contains("swap"))
  cat(paste0(nrow(SNP_df), "rows processed, \n", 
               sum(!iSame_het, na.rm = T), 
               "calls marked as bad_alleles due to allele NAs or mismatches/n"))

  return(SNP_df)
}

ImportBed <- function(bed_file) {
  # bed_file <- repeats_bed
  bed_df <- read.table(bed_file, sep = "\t", header = F)
  colnames(bed_df) <- c("CHROM", "start", "end", "start_POSi", "end_POSi")[1:ncol(bed_df)]
  i_zb_cols <- grepl("start", colnames(bed_df))
  bed_df[, i_zb_cols] <- bed_df[, i_zb_cols] + 1
  if(any(is.numeric(bed_df$CHROM))) {
    bed_df$CHROM <- factor(str_pad(bed_df$CHROM, 2, pad = "0"))
  }
  return(bed_df)
}

MarkRepeats <- function(vcf_in, repeats_in) {
  # vcf_in <- SNPs_merge_finalGT
  # repeats_in$start_POSi <- ConvertPosIndicies(repeats_in, pos_col = "start_POS")
  # repeats_in$end_POSi <- ConvertPosIndicies(repeats_in, pos_col = "end_POS")
  
  rpt <- rep(F, nrow(vcf_in))
  posi <- vcf_in$POSi
  # mrk <- function(x) {
  #   rpt_start <- x$start_POSi[i]
  #   rpt_end <- x$end_POSi[i]
  #   i_rpt <- posi >= rpt_start & posi <= rpt_end
  #   rpt[i_rpt] <- T
  # }

  for(i in 1:nrow(repeats_in)) {
    # i <- 10
    rpt_start <- repeats_in$start_POSi[i]
    rpt_end <- repeats_in$end_POSi[i]
    i_rpt <- posi >= rpt_start & posi <= rpt_end
    rpt[i_rpt] <- T
    if(i %% 1000 == 0) {
      print(paste0(round(i/nrow(repeats_in)*100, 2), "%"))
    }
  }
  return(rpt)
}

GenotypeFromGQ <- function(all_alleleMerge, baseThrsh = 50, 
                           P1_ID = "P1call", P2_ID = "P2call",
                           naThrsh = 100, diffThrsh = 30, 
                           include_DP = T, check_alleles = T) {
  # baseThrsh   baseline threshold - at least one call GQ must exceed this value
  # diffThrsh   difference threshold - if calls do not agree, the higher GQ 
  #             must exceed the lower GQ by this value to be called
  #---------------------------------------------------------------------------#
  # For all sites, the P1-ref-call and P2-ref-call alleles must agree
  # For Ref-hom calls, an Alt allele is not reported, thus only one allele
  # can be compared at these sites. It is also possible to have mutually
  # exclusive allele data, in which both calls are Ref-hom, and thus each
  # provides the allele missing in the partner call. These are a small portion
  # of the data and are eliminated.
  # Split data into sites that agree between references and those that do not
  # Cordant sites do not require a choice in genotype, but should report a
  # level of confidence in terms of GQ value. And so, the higher GQ 
  # value is reported. 
  # Discordant sites must be assigned a genotype from two sets of calls.
  # The GT is awarded to the Ref call with the higher GQ value, so long as it
  # exceeds the lower GQ by diffThrsh. Sites for which the GQ values are too 
  # close are discarded, as little confidence can be rescued with other
  # quality measures.
  # For all sites, the GQ value of the non NA allele must exceed baseThrsh
  # by the naAdjThrsh differential. 
  #---------------------------------------------------------------------------#
  
  # all_alleleMerge <- SNPs_merge %>% filter(ID == "L000")
  random_round <- function(x, seed = 123, tol = 0.1, .digits = 0) { 
    set.seed(seed) 
    round(jitter(x, amount = tol), digits = .digits)
  }
  # Ensure all threshold variables are numeric
  baseThrsh <- as.numeric(baseThrsh)
  naThrsh <- as.numeric(naThrsh) 
  diffThrsh <- as.numeric(diffThrsh)
  
  colnames(all_alleleMerge) <- gsub(P1_ID, "P1call", colnames(all_alleleMerge))
  colnames(all_alleleMerge) <- gsub(P2_ID, "P2call", colnames(all_alleleMerge))
  
  # Include the no-call symbol "./." in the genotype levels in both callsets
  if(any(levels(all_alleleMerge$GT_P1call) == "./.")) {
    GTlvls <- levels(all_alleleMerge$GT_P1call)
  } else {
    GTlvls <- c(levels(all_alleleMerge$GT_P1call), "./.")
  }
  
  # We require the genotype factors to be the same between callsets
  all_alleleMerge$GT_P1call <- factor(all_alleleMerge$GT_P1call, levels = GTlvls)
  all_alleleMerge$GT_P2call <- factor(all_alleleMerge$GT_P2call, levels = GTlvls)
  
  # Index NA rows
  iNA_P1 <- is.na(all_alleleMerge$GQ_P1call) | is.na(all_alleleMerge$GT_P1call)
  iNA_P2 <- is.na(all_alleleMerge$GQ_P2call) | is.na(all_alleleMerge$GT_P2call)
  
  # Convert NA to no-call
  all_alleleMerge$GT_P1call[iNA_P1] <- "./."
  all_alleleMerge$GT_P2call[iNA_P2] <- "./."
  
  # Index NAs for GQ and assign score of 0
  all_alleleMerge$GQ_P1call[iNA_P1] <- 0
  all_alleleMerge$GQ_P2call[iNA_P2] <- 0
  
  
  # All sites are initially assigned missing data values and filled in with 
  # data that passes the thresholds set
  all_alleleMerge$GT <- "./."
  all_alleleMerge$GQ <- 0
  if(include_DP) {
    all_alleleMerge$Ref_DP_final <- 0
    all_alleleMerge$Alt_DP_final <- 0
  }
  # Index same genotypes between callsets
  iNC <- all_alleleMerge$GT_P1call == "./." & all_alleleMerge$GT_P2call == "./."
  iSameGT <- all_alleleMerge$GT_P1call == all_alleleMerge$GT_P2call & !iNC
  
  # Index sites with the same call and for which the P1-ref-call has a higher GQ value 
  iP1overP2 <- iSameGT & all_alleleMerge$GQ_P1call >= all_alleleMerge$GQ_P2call
  
  # Likewise with the P2-ref-call set
  iP2overP1 <- iSameGT & all_alleleMerge$GQ_P2call > all_alleleMerge$GQ_P1call
  
  # Assign values by indicies from above
  all_alleleMerge$GT[iP1overP2] <- as.character(all_alleleMerge$GT_P1call[iP1overP2])
  all_alleleMerge$GT[iP2overP1] <- as.character(all_alleleMerge$GT_P2call[iP2overP1])
  all_alleleMerge$GQ[iP1overP2] <- all_alleleMerge$GQ_P1call[iP1overP2]
  all_alleleMerge$GQ[iP2overP1] <- all_alleleMerge$GQ_P2call[iP2overP1]
  
  if(include_DP) {
    all_alleleMerge$Ref_DP_final[iSameGT] <- random_round(rowSums(all_alleleMerge[iSameGT, c("Ref_DP_P1call", "Ref_DP_P2call")])/2)
    all_alleleMerge$Alt_DP_final[iSameGT] <- random_round(rowSums(all_alleleMerge[iSameGT, c("Alt_DP_P1call", "Alt_DP_P2call")])/2)
  }
  
  # Assign genotypes for discordant sites. The genotype for each site is assigned
  # as either the P1-ref- or P2-ref-derived call, depending on which has a higher
  # GQ value.
  # Index sites for which the P1 call will be used with the logic:
  # genotypes do not agree -and- 
  # # the P1 call is not NA -and- 
  # # # either P2 call is not NA and P1 GQ exceeds P2 by the diffThrsh -or-
  # # # P2 is NA and P1 GQ exceeds the naThrsh
  iP1call <- !iSameGT & 
    !iNA_P1 & 
    ((!iNA_P2 & 
        all_alleleMerge$GQ_P1call > all_alleleMerge$GQ_P2call + diffThrsh) | 
       (iNA_P2 & 
          all_alleleMerge$GQ_P1call >= naThrsh))
  
  iP2call <- !iSameGT & 
    !iNA_P2 &
    ((!iNA_P1 & 
        all_alleleMerge$GQ_P2call > all_alleleMerge$GQ_P1call + diffThrsh) |
       (iNA_P1 & 
          all_alleleMerge$GQ_P2call >= naThrsh))
  
  all_alleleMerge$GT[iP1call] <- as.character(all_alleleMerge$GT_P1call[iP1call])
  all_alleleMerge$GQ[iP1call] <- all_alleleMerge$GQ_P1call[iP1call]
  all_alleleMerge$GT[iP2call] <- as.character(all_alleleMerge$GT_P2call[iP2call])
  all_alleleMerge$GQ[iP2call] <- all_alleleMerge$GQ_P2call[iP2call]
  
  if(include_DP) {
    all_alleleMerge$Ref_DP_final[iP1call] <- all_alleleMerge$Ref_DP_P1call[iP1call]
    all_alleleMerge$Alt_DP_final[iP1call] <- all_alleleMerge$Alt_DP_P1call[iP1call]
    all_alleleMerge$Ref_DP_final[iP2call] <- all_alleleMerge$Ref_DP_P2call[iP2call]
    all_alleleMerge$Alt_DP_final[iP2call] <- all_alleleMerge$Alt_DP_P2call[iP2call]
  }
  
  if(check_alleles) {
  # Revert sites for which alleles do not agree within the set that have received a final call
    all_alleleMerge$REF_P1call[is.na(all_alleleMerge$REF_P1call)] <- "."
    all_alleleMerge$ALT_P1call[is.na(all_alleleMerge$ALT_P1call)] <- "."
    all_alleleMerge$REF_P2call[is.na(all_alleleMerge$REF_P2call)] <- "."
    all_alleleMerge$ALT_P2call[is.na(all_alleleMerge$ALT_P2call)] <- "."
    
    iDiffAllele <-  all_alleleMerge$REF_P1call != "." & 
      all_alleleMerge$REF_P2call != "." &
      !(all_alleleMerge$REF_P1call == all_alleleMerge$REF_P2call &
          all_alleleMerge$ALT_P1call == all_alleleMerge$ALT_P2call)
    # all_alleleMerge[iDiffAllele, ]
    
    all_alleleMerge$GT[iDiffAllele] <- "./."
    all_alleleMerge$GQ[iDiffAllele] <- 0
    if(include_DP) {
      all_alleleMerge$Ref_DP_final[iDiffAllele] <- 0
      all_alleleMerge$Alt_DP_final[iDiffAllele] <- 0
    }
  }
  # Revert sites not exceeding baseThrsh
  iP1underBase <- all_alleleMerge$GQ_P1call < baseThrsh
  iP1underBase[is.na(iP1underBase)] <- T
  iP2underBase <- all_alleleMerge$GQ_P2call < baseThrsh
  iP2underBase[is.na(iP2underBase)] <- T
  iUnderBase <- iP1underBase & iP2underBase
  
  all_alleleMerge$GT[iUnderBase] <- "./."
  all_alleleMerge$GQ[iUnderBase] <- 0
  if(include_DP) {
    all_alleleMerge$Ref_DP_final[iUnderBase] <- 0
    all_alleleMerge$Alt_DP_final[iUnderBase] <- 0
  }
  all_alleleMerge$GT <- factor(all_alleleMerge$GT, levels = GTlvls)
  colnames(all_alleleMerge) <- gsub("P1call", P1_ID, colnames(all_alleleMerge))
  colnames(all_alleleMerge) <- gsub("P2call", P2_ID, colnames(all_alleleMerge))
  return(all_alleleMerge)
}

AddDPinfo <- function(SNPs_finalGT_df) {
  # Takes in a SNP dataframe and appends Sum_DP_final and f_Alt
  # columns, then returns the appended dataframe
  if(!any((colnames(SNPs_finalGT_df) == "Sum_DP_final"))) {
    SNPs_finalGT_df$Sum_DP_final <- SNPs_finalGT_df$Ref_DP_final + SNPs_finalGT_df$Alt_DP_final
  }

  if(!any(colnames(SNPs_finalGT_df) == "f_Alt")) {
    SNPs_finalGT_df$f_Alt <- SNPs_finalGT_df$Alt_DP_final/SNPs_finalGT_df$Sum_DP_final
  }

  return(SNPs_finalGT_df)
}


CalcGenoStats <- function(df_in, group = "all", anc_id = NA) {
  # df_in <- SNPs_merge_finalGT
  
  n_clones <- df_in %>% select(ID) %>% distinct(ID) %>% nrow()
  
  allGT_wide <- df_in %>% 
    ungroup() %>% 
    distinct(ID, POSi, .keep_all = T) %>%
    select(CHROM, POS, POSi, ID, GT) %>% 
    pivot_wider(names_from = ID, values_from = GT)
  
  if(group %in% c("anc", "all")) {
    anc_wide <- allGT_wide %>% select(CHROM, POS, POSi, any_of(anc_id))
    
    n_anc <- ncol(anc_wide) - 3
    n_het <- apply(anc_wide[, -c(1:3)], 1, 
                   function(x) sum(x %in% c("0/1", "het"), na.rm = T))
    n_Ref <- apply(anc_wide[, -c(1:3)], 1, 
                   function(x) sum(x %in% c("0/0", "Ref_hom", "BY"), na.rm = T))
    n_Alt <- apply(anc_wide[, -c(1:3)], 1, 
                   function(x) sum(x %in% c("1/1", "Alt_hom", "RM"), na.rm = T))
    n_Hom <- n_Ref + n_Alt
    n_na <- apply(anc_wide[, -c(1:3)], 1, 
                  function(x) sum(is.na(x)) + sum(x == "./.", na.rm = T))
    anc_wide_vals <- cbind(anc_wide[, 1:3], nRef_anc = n_Ref, 
                           nHet_anc = n_het, nAlt_anc = n_Alt, nHom_anc = n_Hom, nNA_anc = n_na)
    
    anc_wide_vals$fRef_anc <- anc_wide_vals$nRef_anc/(n_anc - anc_wide_vals$nNA_anc)
    anc_wide_vals$fHet_anc <- anc_wide_vals$nHet_anc/(n_anc - anc_wide_vals$nNA_anc)
    anc_wide_vals$fAlt_anc <- anc_wide_vals$nAlt_anc/(n_anc - anc_wide_vals$nNA_anc)
    
    if(group == "anc") {
      return(anc_wide_vals)
      stop("ancestor table returned")
    }
  }
  if(group %in% c("evo", "all")) {
    if(group == "all") {
      evoHet_wide <- allGT_wide %>% select(CHROM, POS, POSi, !any_of(anc_id))
    } else {
      evoHet_wide <- allGT_wide
    }
    
    n_evo <- ncol(evoHet_wide) - 3
    n_het <- apply(evoHet_wide[, -c(1:3)], 1, 
                   function(x) sum(x %in% c("0/1", "het"), na.rm = T))
    n_Ref <- apply(evoHet_wide[, -c(1:3)], 1, 
                   function(x) sum(x %in% c("0/0", "Ref_hom", "BY"), na.rm = T))
    n_Alt <- apply(evoHet_wide[, -c(1:3)], 1, 
                   function(x) sum(x %in% c("1/1", "Alt_hom", "RM"), na.rm = T))
    n_na <- apply(evoHet_wide[, -c(1:3)], 1, 
                  function(x) sum(is.na(x)) + sum(x == "./.", na.rm = T))
    evo_wide_vals <- cbind(evoHet_wide[, 1:3], nRef_evo = n_Ref, 
                           nHet_evo = n_het, nAlt_evo = n_Alt, nNA_evo = n_na)
    evo_wide_vals$fRef_evo <- evo_wide_vals$nRef_evo/(n_evo - evo_wide_vals$nNA_evo)
    evo_wide_vals$fHet_evo <- evo_wide_vals$nHet_evo/(n_evo - evo_wide_vals$nNA_evo)
    evo_wide_vals$fAlt_evo <- evo_wide_vals$nAlt_evo/(n_evo - evo_wide_vals$nNA_evo)
    
    if(group == "evo") {
      evo_wide_vals <- evo_wide_vals %>% 
        # rename_with(.fn = ~paste0(.x, "_evo"),
        #             .cols = !c(CHROM, POS, POSi)) %>%
        mutate(nTot_evo = n_clones - nNA_evo,
               convBias_evo = nRef_evo/(nRef_evo + nAlt_evo),
               nConv_evo = nRef_evo + nAlt_evo)
      
      return(evo_wide_vals)
      stop("end-point table returned")
    }
  }
  
  if(group == "all") {
    all_wide_vals <- merge(anc_wide_vals, evo_wide_vals, 
                           by = c("CHROM", "POS", "POSi"), 
                           all = T) %>% 
      arrange(POSi) %>%
      mutate(nTot_evo = n_clones - nNA_evo,
             convBias_evo = nRef_evo/(nRef_evo + nAlt_evo),
             nConv_evo = nRef_evo + nAlt_evo)
    
    return(all_wide_vals)
    stop("all clones table returned")
  }
}

BB_alphaFromVar <- function(.v, .n) {
    # https://en.wikipedia.org/wiki/Beta-binomial_distribution
    # Variance equation, solve for alpha when alpha = beta
    .a <- (((.n - 1)/(4 * .v / .n - 1)) - 1) / 2
    return(.a)
}

BB_alphaFromDP <- function(.n, .b) {
  # Combines formula of BB_alphaFromVar with linear model
  # of variance vs depth under binomial + x^2 error term
  .a = (0.25/.b - 0.25/(.b * .n) - 1)/2
  return(.a)
}

qbetabinom <- function(.x, .prob = NULL, .size, .theta = NULL, .alpha = NULL, .beta = NULL, .p) {
  # .x <- 20
  # .prob <- 0.5
  # .size <- .x
  # .theta <- .x * 4
  # .alpha <- 58; .beta <- .alpha
  # .p <- 0.05
  if(is.null(.theta)) {
    if(is.null(.beta)) {
      .beta <- .alpha
    }
    dbb_df <- data.frame(x = 0:.x,
                         p = cumsum(dbetabinom(x = 0:.x, size = .size, shape1 = .alpha, shape2 = .beta)))
  } else {
    dbb_df <- data.frame(x = 0:.x,
                         p = cumsum(dbetabinom(x = 0:.x, prob = .prob, size = .size, theta = .theta)))
  }
  # dbb_df <- data.frame(x = 0:.x, 
  #                      p = cumsum(dbetabinom(x = 0:.x, size = .x,  shape1 = sqrt(.x), shape2 = sqrt(.x))))
  # 
  if(.p < 0.5) {
    dbb_df_lt <- dbb_df %>% filter(x < .x / 2)
    quant <- max(dbb_df_lt$x[dbb_df_lt$p < .p])
    quant <- ifelse(quant < 0 | is.infinite(quant), 0, quant)
  } else {
    dbb_df_rt <- dbb_df %>% filter(x > .x / 2)
    quant <- min(dbb_df_rt$x[dbb_df_rt$p > .p]) + 1
    quant <- ifelse(quant > .x | is.infinite(quant), .x, quant)
  }
  return(quant)
}

BB_rejectHets <- function(GT_AltDP_SumDP, error_beta = NULL, BB_alpha = NULL, FE_rate, min_DP = 6) {
  # Final output will be a logical column of rejected het calls under the
  # BB model and a given .beta, which is the error parameter of Sum_DP^2 
  # Column is formed on input df, then hets are indexed and
  # a second df is produced for het filtering
  # GT_AltDP_SumDP <- SNPs_merge_finalGT %>% filter(GT == "0/1", Sum_DP_final >= 4) %>% slice_sample(n = 10000) %>% select(GT, Alt_DP_final, Sum_DP_final)
  # error_beta <- 0.00291
  colnames(GT_AltDP_SumDP) <- c("GT", "Alt_DP", "Sum_DP")
  GT_AltDP_SumDP$BB_reject <- F
  i_het <- GT_AltDP_SumDP$GT == "0/1" & GT_AltDP_SumDP$Sum_DP >= min_DP
  i_het_lo_DP <- GT_AltDP_SumDP$GT == "0/1" & GT_AltDP_SumDP$Sum_DP < min_DP
  BB_het_df <- GT_AltDP_SumDP[i_het, ]
  if((!is.null(error_beta) & !is.null(BB_alpha)) | 
     (is.null(error_beta) & is.null(BB_alpha))) {
    cat("Must input one parameter, error_beta or BB_alpha")
  } else if(!is.null(error_beta)) {
    # define a function that outputs a BB alpha parameter
    # for each Sum_DP and fixed .beta and calculate from het Sum_DPs
    BB_alphaFromDP <- function(dp, .beta) {
      .alpha = (0.25/.beta - 0.25/(.beta * dp) - 1)/2
      return(.alpha)
    }
    BB_het_df$BB_alpha <- BB_alphaFromDP(BB_het_df$Sum_DP, 
                                            .beta = error_beta)
  } else if(!is.null(BB_alpha)) {
    BB_het_df$BB_alpha <- BB_alpha
  }
  
  # Mark rejected het calls by calculating the BB thresholds for each row
  # using mapply to combine each Sum_DP and alpha. FE_rate is the tolerance
  # in rejecting true heterozygous calls

  BB_het_df <- BB_het_df %>%
              group_by(Sum_DP) %>%
              mutate(BB_thresh_hi = qbetabinom(mean(Sum_DP), 0.5, mean(Sum_DP), .alpha = mean(BB_alpha), .p = 1 - FE_rate/2),
                      BB_thresh_lo = qbetabinom(mean(Sum_DP), 0.5, mean(Sum_DP), .alpha = mean(BB_alpha), .p = FE_rate/2)) %>%
              ungroup() %>%
              mutate(BB_reject = Alt_DP > BB_thresh_hi | Alt_DP < BB_thresh_lo)
  # Overwrite the BB_reject column of the original df at het calls
  # with that determined by the model
  GT_AltDP_SumDP$BB_reject[i_het] <- BB_het_df$BB_reject
  GT_AltDP_SumDP$BB_reject[i_het_lo_DP] <- T
  return(GT_AltDP_SumDP$BB_reject)
}

lambdaFromFractK <- function(emp_dist, k, max_lambda = 1.5) {
  # emp_dist <- hom_sites_AD %>% filter(Sum_DP == 236) %>% pull(Alt_DP)
  max_emp <- max(emp_dist)
  n_emp <- length(emp_dist)
  emp_hist <- hist(emp_dist, plot = F, breaks = (0:(max_emp + 1) - 0.5))
  # emp_hist$mids are the Alt DP bins
  if(k > max(emp_hist$mids)) {
    return(NA)
  } else {
    # emp_fraction is the fraction of the dist in each Alt DP bin
    emp_fraction <- emp_hist$counts/n_emp
    f_K <- emp_fraction[emp_hist$mids == k]
    if(f_K == 0) {
      return(NA)
    } else {
      if(k == 0) {
        return(-log(f_K))
      } else {
        # We want to estimate the Lambda of a Poisson that would be
        # expected from a given faction of the data in Alt DP bin k 
        rootLambda = function(x) {
          f_K - x^k * exp(-x) / factorial(k)
        }
        # 
        if((k == 1 & f_K > 0.3678794) | (k == 2 & f_K > 0.2706706)) {
          return(NA)
        } else {
          return(uniroot(rootLambda, interval = c(0, max_lambda))$root)
        }
      }
    }
  }
}

LambdaFromsetKs <- function(emp_dist, k = 0:2) {
  max_emp <- max(emp_dist)
  n_emp <- length(emp_dist)
  emp_hist <- hist(emp_dist, plot = F, breaks = (0:(max_emp + 1) - 0.5))
  mean_emp <- mean(emp_dist)
 
 Ls <- sapply(k, function(x) 
      lambdaFromFractK(emp_dist = emp_dist, k = x, max_lambda = ifelse(mean_emp < 1, 1, mean_emp)))
    return(mean(Ls))
    
}

LambdaLMthrsh <- function(dp, lin_param, exp_param, .alpha = NA, alt_dp = NA) {
    l <- exp(lin_param) * dp ^ exp_param
    if(is.na(.alpha) & !is.na(alt_dp)) {
      .alpha <- 1 - sum(dpois(0:alt_dp, lambda = l))
      l_thrsh <- qpois(p = .alpha, lambda = l)
      return(l_thrsh)
  } else if(!is.na(.alpha) & is.na(alt_dp)) {
    l_thrsh <- qpois(p = .alpha, lambda = l)
    return(l_thrsh)
  } else {
    cat("Provide .alpha OR alt_DP")
  }
}