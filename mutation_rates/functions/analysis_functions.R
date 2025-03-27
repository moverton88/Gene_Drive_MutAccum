
mean.fun <- function(d, idx) {
  mean((d[idx]), trim = 0, na.rm = T)
}

sum.fun <- function(d, idx) {
  sum((d[idx]), trim = 0, na.rm = T)
}

proportion_z_score <- function(c_table, grp) {
  # c_table <- LOH_type_xTx[1:3]
  # grp <- "Tx_ID"
  c_table[, grp] <- factor(c_table[, grp])
  grp_levels <- levels(c_table[, grp])
  combos <- expand.grid(as.numeric(c_table[, grp]), 
                        as.numeric(c_table[, grp])) %>% 
    filter(Var1 < Var2)
  c_1_range <- unique(combos$Var1) %>% sort()
  c_2_range <- unique(combos$Var2) %>% sort()
  Vars <- colnames(c_table)[colnames(c_table) != grp]
  c_table$sums <- c_table[,Vars[1]] + c_table[, Vars[2]]
  
  c_stats <- data.frame(x_1 = c_table[combos$Var1, Vars[1]], 
                        x_2 = c_table[combos$Var2, Vars[1]], 
                        n_1 = c_table$sums[combos$Var1], 
                        n_2 = c_table$sums[combos$Var2])
  c_stats <- c_stats %>% mutate(p_1 = x_1/n_1, p_2 = x_2/n_2)
  combos <- cbind(combos, c_stats)
  combos$z <- 0
  for(c in 1:nrow(combos)) {
    # c <- 1
    c_1 <-  combos$Var1[c]
    c_2 <-  combos$Var2[c]
    p_a <- (combos$x_1[c] + c_stats$x_2[c])/(c_stats$n_1[c] + c_stats$n_2[c])
    combos$z[c] <- (combos$p_1[c] - combos$p_2[c]) /
      sqrt(p_a*(1-p_a)*(1/combos$n_1[c] + 1/combos$n_2[c]))
  }
  combos$p <- pnorm(q = -abs(combos$z), lower.tail = T) * 2
  combos$Var1 <- factor(combos$Var1, labels = grp_levels[c_1_range])
  combos$Var2 <- factor(combos$Var2 - 1, labels = grp_levels[c_2_range])
  combos_out <- combos %>% select(Var2, Var1, z, p)
  colnames(combos_out) <- c(paste0(grp, "_null"), paste0(grp, "_test"), "Z-score", "P-value")
  return(combos_out)
}

format_sci_10 <- function(l) {
  t_exp <- ifelse(l == 0, 0, 
                  parse(text = gsub("e", " %*% 10^", format(x = l, digits = 2, scientific = T))))
}

factorIntegers <- function(df_variable) {
  # Find length of largest element and pad fewer digit integers
  # with "0"
  max_pad <- nchar(max(df_variable))
  df_variable <- factor(str_pad(df_variable, width = max_pad, pad = "0"))
}

operate_by_factor_match <- function(op_pair, trgt_pair, .fun) {
  # op_pair <- chrom_arms[, c("CHROM", "arm_2_cent")]
  # trgt_pair <- POSi_data_in[!i_cent_neg, c("CHROM", "est_mid_POS")]
  op_pair <- as.data.frame(op_pair)
  trgt_pair <- as.data.frame(trgt_pair)
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
      return(trgt_out)
    } else {
      print("factor levels do not match")
    }
  } else {
    print("factor levels are different lengths")
  }
}

Recode_Tx <- function(Tx_col) {
  # Input the DF variable with one-character 
  # Tx identifier and output a variable with
  # full treatment names
  require(DescTools)
  Tx_name_out <- Recode(Tx_col, "WT" = "N", "Cas9" = "H", "Drive" = "F")
  Tx_name_out <- factor(Tx_name_out, levels = c("WT", "Cas9", "Drive"))
  return(Tx_name_out)
}

Recode_GT <- function(GT_in, RMcall_rev = FALSE) {
  require(DescTools)
  if(RMcall_rev == T) {
    GT_tmp <- Recode(GT_in, "RM_hom" = "0/0", "het" = "0/1", "BY_hom" = "1/1")
    GT_out <- factor(GT_tmp, levels = c("BY_hom", "het", "RM_hom"))
  } else if(RMcall_rev == F) {
    GT_tmp <- Recode(GT_in, "BY_hom" = "0/0", "het" = "0/1", "RM_hom" = "1/1")
    GT_out <- factor(GT_tmp, levels = c("BY_hom", "het", "RM_hom"))
  } else {
    print("RMcall_rev must be a logical")
  }
  return(GT_out)
}

factorTxNames <- function(TxNames_col) {
  TxNames_out <- factor(TxNames_col, levels = c("WT", "Cas9", "Drive"))
  return(TxNames_out)
}

Recode_Tx_ID <- function(Tx_col, tx_type = "Tx", Tx_ID_lvls = Tx_ID_levels) {
  # Input the DF column with old one-character 
  # Tx identifier and output a factored variable with
  # new treatment identifier
  require(DescTools)
  if(!is.factor(Tx_col)) {
    Tx_col <- factor(Tx_col)
  }
  if(tx_type == "Tx") {
    Tx_ID_out <- Recode(Tx_col, "W" = "N", "C" = "H", "D" = "F")
  }
   else if(tx_type == "Tx_name"){
    Tx_ID_out <- Recode(Tx_col, "W" = "WT", "C" = "Cas9", "D" = "Drive")
   } else if(tx_type == "Tx_ID"){
     Tx_ID_out <- Tx_col
   } else {
    print("Treatment levels not recognized")
  }
  Tx_ID_out <- factor(Tx_ID_out, levels = Tx_ID_lvls)
  return(Tx_ID_out)
}

CategoriesFromID <- function(df_in) {
  # requires variable "ID" in dataframe with the format T_Lxx
  # where T is the Tx, T_L is the Line, and xx is the Rep
  df_in$Tx <- factor(substr(df_in$ID, 1, 1))
  df_in$Tx_ID <- Recode_Tx_ID(df_in$Tx, tx_type = "Tx")
  df_in$Line <- factor(substr(df_in$ID, 1, 3))
  df_in$Rep <- factor(substr(df_in$ID, 4, 5))
  df_in$Tx_name <- Recode_Tx(df_in$Tx)
  return(df_in)
}


chainToDF <- function(chainFile, chrom_bound_1 = chrom_bound_BY, chrom_bound_2 = chrom_bound_RM)  {
  chainTable <- read_lines(chainFile)
  chainTable <- gsub(pattern = "\"", replacement = "", x = chainTable)
  chainTable <- chainTable[!chainTable == ""]
  chrIdx <- which(chainTable %in% chrom_bound_BY$CHROM)
  endIdx <- length(chainTable)
  valueIdx <- data.frame(from = chrIdx + 1, to = c(chrIdx[2:16] - 1, endIdx))
  chrNames <- chainTable[chrIdx]
  chainDf <- data.frame(NULL)
  for (ch in 1:length(chrNames)) {
    # ch = 3
    # subName <- chainTable[chrIdx[ch]]
    subChain <- chainTable[valueIdx$from[ch]:valueIdx$to[ch]]
    subDF <- colsplit(subChain, " ", names = c("block", "BYdiff", "RMdiff"))
    subDF[is.na(subDF)] <- 0
    subDF$csBlock <- cumsum(subDF$block)
    subDF$csBYdiff <- cumsum(subDF$RMdiff) - cumsum(subDF$BYdiff)
    subDF$csRMdiff <- cumsum(subDF$BYdiff) - cumsum(subDF$RMdiff)
    subDF$BY_POS <- subDF$csBlock + lag(cumsum(subDF$BYdiff)) + 1
    subDF$RM_POS <- subDF$csBlock + lag(cumsum(subDF$RMdiff)) + 1
    subDF[1, c("BY_POS", "RM_POS")] <- subDF$csBlock[1] + 1
    subDF$BY_POSi <- subDF$BY_POS + chrom_bound_1$Start[ch]
    subDF$RM_POSi <- subDF$RM_POS + chrom_bound_2$Start[ch]
    subDF$CHROM <- chrNames[ch]
    chainDf <- rbind(chainDf, subDF)
    chainDf$CHROM <- factor(chainDf$CHROM)
  }
  return(chainDf)
}

BoundsFromLengths <- function(lengths) {
  # lengths <- chrom_lengths_P1$length
  chr_start <- c(1, cumsum(lengths[1:(length(lengths) - 1)]) + 1)
  chr_end <- cumsum(lengths)

  chr_bound <- data.frame(CHROM = factor(str_pad(1:length(lengths), width = 2, pad = "0")), 
                          Start = as.numeric(chr_start), 
                          End = as.numeric(chr_end))
  return(chr_bound)
}

GetChromDF <- function(chrom_lengths) {
  # chrom_lengths <- chrom_lengths_Ref
  n_chroms <- nrow(chrom_lengths)
  roman_chr <- factor(as.character(as.roman(1:n_chroms)), levels = as.character(as.roman(1:n_chroms)))
  chrom_IDs <- data.frame(CHROM = factor(str_pad(1:n_chroms, 2, pad = "0")), rom_CHROM = roman_chr)
  chrom_lengths <- merge(chrom_IDs, chrom_lengths, by = "CHROM")
  
  chrom_indcs <- cumsum(chrom_lengths$length)
  
  chrom_bounds <- BoundsFromLengths(chrom_lengths$length)
  chrom_bounds <- cbind(chrom_lengths, chrom_bounds[c("Start", "End")])
  
  return(chrom_bounds)
}


ChainToDFlist <- function(chain_file, remove_mito = F)  {
  ###########################################################################
  # A chain file is composed of a series of sections, one for each contig or 
  # chromosome. Within each section, there are three columns:
  # 1) block length, 
  # 2) difference between end of this block and start of next block in the 
  # reference sequence, 
  # 3) difference between end of this bloc and beginning of next query block
  # This function reads in the chain file and converts it to a data.frame 
  # capable of converting the position coordinates of the query vcf back to 
  # that of the coordinates of the Type reference sequence.
  # It also creates data.frames of the chromosome lengths
  
  # chain_file <- P2_chainFile
  chain_table <- read_lines(chain_file) # read in chain file
  chain_table <- gsub(pattern = "\"", replacement = "", x = chain_table) # remove spurious "\" characters
  chain_table <- chain_table[!chain_table == ""] # remove new lines at end of each section
  chr_idx <- grep("chain", chain_table) # index section headers
  if(remove_mito) {
    # remove mitochondrial sequence if present and is present as the final section of the chain file
    mito_idx <- chr_idx[length(chr_idx)]
    chain_table <- chain_table[1:(mito_idx - 1)]
    chr_idx <- chr_idx[1:(length(chr_idx) - 1)]
  } 
  end_idx <- length(chain_table)
  value_idx <- data.frame(from = chr_idx + 1, to = c(chr_idx[2:length(chr_idx)] - 1, end_idx)) # dataframe of section blocks
  chr_headers <- chain_table[chr_idx] # get section headers
  chr_info <- str_split(chr_headers, " ", simplify = T) %>% as.data.frame() # turn headers into dataframe
  col_names <- c("CHROM", "length", "strand", "start", "end") # rename dataframe colnames
  colnames(chr_info) <- c("null", "score", paste0("Ref_", col_names), paste0("Qry_", col_names), "ID")
  chrom_names <- str_pad(1:nrow(chr_info), width = 2, side = "left", pad = "0") # get chromosome names
  
  Ref_chrom_starts <- cumsum(c(1, chr_info$Ref_length[1:(nrow(chr_info) - 1)]))
  Qry_chrom_starts <- cumsum(c(1, chr_info$Qry_length[1:(nrow(chr_info) - 1)]))
  
  chain_df <- data.frame(NULL)
  for (ch in 1:nrow(chr_info)) {
    # ch = 1
    # subName <- chain_table[chr_idx[ch]]
    sub_chain <- chain_table[value_idx$from[ch]:value_idx$to[ch]]
    sub_df <- colsplit(sub_chain, " ", names = c("block", "Ref_diff", "Qry_diff"))
    sub_df[is.na(sub_df)] <- 0
    sub_df$cs_block <- cumsum(sub_df$block)
    sub_df$cs_Ref_diff <- cumsum(sub_df$Qry_diff) - cumsum(sub_df$Ref_diff)
    sub_df$cs_Qry_diff <- cumsum(sub_df$Ref_diff) - cumsum(sub_df$Qry_diff)
    sub_df$Ref_POS <- sub_df$cs_block + cumsum(sub_df$Ref_diff) + 1 # cs_block is 0 indexed, so add 1 for position
    sub_df$Qry_POS <- sub_df$cs_block + cumsum(sub_df$Qry_diff) + 1
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

RMxBY_liftover <- function(vcf_df, chain_df = BYtoRMchainDf, 
                           lift_from = "RM", lift_to = "BY",
                           chrom_from = chrom_bound_RM,
                           chrom_to = chrom_bound_BY) {
  require(dplyr)
  # vcf_df <- all_fltrd_RM %>% filter(ID == "N_A00")
  if(!(lift_to %in% c("BY", "RM"))) {
    print("lift_to must be 'BY' or 'RM'")
  } else {
    posi_col <- paste0(lift_from, "_POSi") 
    csDiff_col <- paste0("cs", lift_from, "diff")
  }
  vcf_df$POS_adj <- vcf_df$POS
  vcf_df$POSi_adj <- vcf_df$POSi
  # chain_df$gDiff <- chain_df[, gDiff_col_1] - chain_df[, gDiff_col_2]
  for(ch in levels(chain_df$CHROM)) {
    # ch = levels(BYtoRMchainDf$CHROM)[3]
    ch_chain <- chain_df %>% filter(CHROM == ch)
    chDiff <- chrom_to$Start[chrom_to$CHROM == ch] - chrom_from$Start[chrom_from$CHROM == ch]
    # iChrom_vcf <- vcf_df$CHROM == ch
    nBlocks <- nrow(ch_chain)
    for(b in 1:(nBlocks - 1)) {
      # b = 1
      b_i <- vcf_df$POSi >= ch_chain[b, posi_col] & vcf_df$POSi < ch_chain[(b + 1), posi_col]
      vcf_df$POS_adj[b_i] <- vcf_df$POS[b_i] + ch_chain[b, csDiff_col]
      vcf_df$POSi_adj[b_i] <- vcf_df$POSi[b_i] + ch_chain[b, csDiff_col] + chDiff
    }
    e_i <- vcf_df$POSi >= ch_chain[nBlocks, posi_col] & vcf_df$POSi <= chrom_from$End[as.numeric(ch)]
    vcf_df$POS_adj[e_i] <- vcf_df$POS[e_i] + ch_chain[nBlocks, csDiff_col]
    lDiff <-  chrom_to$Start[nrow(chrom_to)] - chrom_from$Start[nrow(chrom_to)]
    vcf_df$POSi_adj[e_i] <- vcf_df$POSi[e_i] + ch_chain[nBlocks, csDiff_col] + lDiff
    print(ch)
  }
  vcf_df$POS <- vcf_df$POS_adj
  vcf_df$POSi <- vcf_df$POSi_adj
  vcf_df <- vcf_df %>% dplyr::select(-c("POS_adj", "POSi_adj"))
  return(vcf_df)
}

SwapCrossedAlleles <- function(SNP_df, annotated_POSi, to_match = "BYcall",
                               to_swap = "RMcall") {
  # We expect crossed alleles at sites that differ between the reference sequences.
  # This should be the same set of sites as is present in RMxBY_vcf
  
  # SNP_df <- SNPs_merge
  # annotated_POSi <- RMxBY_comp_SNPs_POSi
  
  # Get columns to match and to swap
  cnm <- colnames(SNP_df)
  
  i_match <- grep(to_match, cnm)
  i_swap <- grep(to_swap, cnm)
  i_ref_alel <- grep("REF", cnm)
  i_alt_alel <- grep("ALT", cnm)
  
  i_alel_ord <- c(intersect(i_match, i_ref_alel), 
                  intersect(i_match, i_alt_alel),
                  intersect(i_swap, i_ref_alel),
                  intersect(i_swap, i_alt_alel))
  
  i_ref_DP <- grep("Ref_DP", cnm)
  i_alt_DP <- grep("Alt_DP", cnm)
  i_DP_ord <- c(intersect(i_match, i_ref_DP), 
                intersect(i_match, i_alt_DP),
                intersect(i_swap, i_ref_DP),
                intersect(i_swap, i_alt_DP))
  
  i_GT <- grep("GT", cnm)
  i_GT_ord <- c(intersect(i_match, i_GT), 
                intersect(i_swap, i_GT))
  
  if(length(c(i_alel_ord, i_DP_ord, i_GT_ord)) != 10) {
    print("data does not contain the correct columns")
  }
  
  # Index sites where the REF allele of one callset matches the ALT allele of the other
  iCross_alleles <- SNP_df[, i_alel_ord[1]] == SNP_df[, i_alel_ord[4]] & 
    SNP_df[, i_alel_ord[2]] == SNP_df[, i_alel_ord[3]]
  
  # Do not include sites with NAs by assigning as FALSE
  iCross_alleles[is.na(iCross_alleles)] <- F
  
  i_RMxBY_merge <- SNP_df$POSi %in% annotated_POSi
  # i_RMxBY_merge[3:4] <- T
  
  # Record POSi for which the alleles are crossed but missing in the RMxBY df or
  # for which they are present in the RMxBY df but not crossed as expected
  i_dscd_source <- (!i_RMxBY_merge & iCross_alleles) | (i_RMxBY_merge & !iCross_alleles)
  
  # The set of sites to swap alleles: Known from references or data and not 
  # missing in the RM callset
  iSwap_alleles <- iCross_alleles & !is.na(SNP_df$REF_RMcall)
  
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
  
  # Index homozygous sites in RM callset
  iAlt_RMcall <- SNP_df[, i_GT_ord[2]] == "0/0"
  iRef_RMcall <- SNP_df[, i_GT_ord[2]] == "1/1"
  
  # Edit original GT column to swap genotypes of crossed alleles
  SNP_df[iSwap_alleles & iAlt_RMcall, i_GT_ord[2]] <- "1/1"
  SNP_df[iSwap_alleles & iRef_RMcall, i_GT_ord[2]] <- "0/0"
  
  # Index matching alleles after swap
  iSame_het <- SNP_df[, i_alel_ord[1]] == SNP_df[, i_alel_ord[3]] & SNP_df[, i_alel_ord[2]] == SNP_df[, i_alel_ord[4]]
  
  # Column of sites which still do not match correctly
  SNP_df$bad_alleles <- !iSame_het | i_dscd_source
  
  SNP_df <- SNP_df %>% select(!contains("swap"))
  print(paste0(nrow(SNP_df), "rows processed, ", 
               sum(!iSame_het), 
               "rows marked as bad_alleles due to persistent allele mismatches"))
  return(SNP_df)
}

ChromosomeCoordinates <- function(df_in, chrom_bound = chrom_bound_BY, POSi_col = "POSi") {
  # df_in <- mean_LOHrates50k
  df_in <- as.data.frame(df_in)
  df_in$CHROM <- factor(df_in$CHROM)
  df_in$dist_cent <- 0
  df_in$dist_term_1 <- 0
  df_in$dist_term_2 <- 0
  for(chr in centrom_df$CHROM){
    # chr <- "02"
    df_in$dist_cent[df_in$CHROM == chr] <-
      df_in[df_in$CHROM == chr, POSi_col] - 
      centrom_df$POSi[centrom_df$CHROM == chr]
    df_in$dist_term_1[df_in$CHROM == chr] <-
      df_in[df_in$CHROM == chr, POSi_col] - 
      chrom_bound$Start[chrom_bound$CHROM == chr]
    df_in$dist_term_2[df_in$CHROM == chr] <-
      df_in[df_in$CHROM == chr, POSi_col] - 
      chrom_bound$End[chrom_bound$CHROM == chr]
  }
  df_in$dist_term <- ifelse(df_in$dist_cent < 0, df_in$dist_term_1, - df_in$dist_term_2)
  df_in$dist_fract <- abs(df_in$dist_cent) / (abs(df_in$dist_cent) + df_in$dist_term)
  df_in$abs_dist_fract <- abs(df_in$dist_cent) / (abs(df_in$dist_cent) + df_in$dist_term)
  return(df_in)
}

ConvertPosIndicies <- function(pos_df, pos_col = "POS", chrom_col = "CHROM", 
                                        index_out = "POSi", chrom_bound = chrom_bound_BY,
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

MarkRepeats <- function(vcf_in, repeats_in = repeats_bed) {
  # vcf_in <- SNPs_merge_finalGT
  repeats_in <- as.data.frame(repeats_in)
  repeats_in$Start_POSi <- ConvertPosIndicies(repeats_in, pos_col = "Start_POS")
  repeats_in$End_POSi <- ConvertPosIndicies(repeats_in, pos_col = "End_POS")
  
  rpt <- rep(F, nrow(vcf_in))
  posi <- vcf_in$POSi
  for(i in 1:nrow(repeats_in)) {
    # i <- 10
    rpt_start <- repeats_in$Start_POSi[i]
    rpt_end <- repeats_in$End_POSi[i]
    i_rpt <- posi >= rpt_start & posi <= rpt_end
    rpt[i_rpt] <- T
    if(i %% 100 == 0) {
      print(paste0(round(i/nrow(repeats_in)*100, 2), "%"))
    }
  }
  return(rpt)
}

site_genotype_stats <- function(df_in, group = "all") {
  # df_in <- SNPs_noAnc
  allGT_wide <- df_in %>% 
    ungroup() %>% 
    distinct(ID, POSi, .keep_all = T) %>%
    select(CHROM, POS, POSi, ID, GT) %>% 
    pivot_wider(names_from = ID, values_from = GT)
  
  if(group %in% c("anc", "all")) {
    anc_wide <- allGT_wide %>% select(CHROM, POS, POSi, contains("00"))
    
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
    anc_wide_vals <- cbind(anc_wide[, 1:3], nRef = n_Ref, 
                           nHet = n_het, nAlt = n_Alt, nHom_anc = n_Hom, nNA = n_na)
    
    anc_wide_vals$fRef <- anc_wide_vals$nRef/(n_anc - anc_wide_vals$nNA)
    anc_wide_vals$fHet <- anc_wide_vals$nHet/(n_anc - anc_wide_vals$nNA)
    anc_wide_vals$fAlt <- anc_wide_vals$nAlt/(n_anc - anc_wide_vals$nNA)
    if(group == "anc") {
      return(anc_wide_vals)
      stop("ancestor table returned")
    }
  }
  if(group %in% c("evo", "all")) {
    evoHet_wide <- allGT_wide %>% select(CHROM, POS, POSi, !contains("00"))
    
    n_evo <- ncol(evoHet_wide) - 3
    n_het <- apply(evoHet_wide[, -c(1:3)], 1, 
                   function(x) sum(x %in% c("0/1", "het"), na.rm = T))
    n_Ref <- apply(evoHet_wide[, -c(1:3)], 1, 
                   function(x) sum(x %in% c("0/0", "Ref_hom", "BY"), na.rm = T))
    n_Alt <- apply(evoHet_wide[, -c(1:3)], 1, 
                   function(x) sum(x %in% c("1/1", "Alt_hom", "RM"), na.rm = T))
    n_na <- apply(evoHet_wide[, -c(1:3)], 1, 
                  function(x) sum(is.na(x)) + sum(x == "./.", na.rm = T))
    evo_wide_vals <- cbind(evoHet_wide[, 1:3], nRef = n_Ref, 
                           nHet = n_het, nAlt = n_Alt, nNA = n_na)
    evo_wide_vals$fRef <- evo_wide_vals$nRef/(n_evo - evo_wide_vals$nNA)
    evo_wide_vals$fHet <- evo_wide_vals$nHet/(n_evo - evo_wide_vals$nNA)
    evo_wide_vals$fAlt <- evo_wide_vals$nAlt/(n_evo - evo_wide_vals$nNA)
    if(group == "evo") {
      return(evo_wide_vals)
      stop("end-point table returned")
    }
  }
  all_wide_vals <- merge(anc_wide_vals, evo_wide_vals, 
                         by = c("CHROM", "POS", "POSi"), 
                         all = T, suffixes = c("_anc", "_evo")) %>% arrange(POSi)
  if(group == "all") {
    return(all_wide_vals)
    stop("all clones table returned")
  }
}

nSites_Calc <- function(all_vcf_in, min_sites = 20000, cvr_plot = F) {
  ## Takes dataframe containing all clones and counts the number of valid sites
  ## in each. It also produces a point plot of these values for each clone.
  require(dplyr)
  require(ggplot2)
  # Count the number of sites in each clone
  ID_siteCounts <- all_vcf_in %>% 
    dplyr::group_by(Tx_name, Line, Rep) %>%
    dplyr::count(ID)
  
  if (cvr_plot == T) {
  # Plot number of valid sites in each sample
  
  cvrgPlot <- ggplot() + 
    geom_jitter(data=subset(ID_siteCounts, Rep != "00"), 
                aes(x=Tx, y=n, group=ID), size=2, height=0, width=0.3, alpha=0.8) + 
    geom_jitter(data=subset(ID_siteCounts, Rep == "00"), 
                aes(x=Tx, y=n, group=ID), color="orange1", size=2,  height=0, width=0.3) + 
    geom_hline(aes(yintercept=min_sites), linetype="dashed") +
    geom_text_repel(data=subset(ID_siteCounts, n < min_sites), 
                    aes(x=Tx, y=n, label=ID, group=ID), point.padding=0.2, force=0.5, segment.colour = NA) +
    xlab("Drive Construct") + ylab("Number of valid sites") 
  
  ggsave(file.path(outIntDir, "cvrgPlot_xTx.png"), 
         plot = cvrgPlot,
         device = "png",
         width = 11, height = 8.5, 
         units = "in",
         dpi = 300)
  
  }
  return(ID_siteCounts)
}

nSites_Filter <- function(all_vcf_in, siteCnts, min_sites = 20000) {
  ## Filters dataframe to only clones with called sites < min_sites
  ## "siteCnts" is the output of the "LOH_cvrgCalc" function
  require(dplyr)
  ID_cut <- siteCnts %>% filter(n < min_sites)
  all_hiCvrg <- all_vcf_in %>% filter(!(ID %in% ID_cut$ID))
  return(all_hiCvrg)
}

filters_failed <- function(df, filters) {
  # df is the vcf df to obtain filter values and provide filter results for
  filter_length <- length(filters)
  filter_split <- (str_split(filters, " "))
  filter_names <- unlist(lapply(filter_split, "[[", 1))
  filter_ops <-  unlist(lapply(filter_split, "[[", 2))
  filter_values <-  as.numeric(unlist(lapply(filter_split, "[[", 3)))
  filter_states <- rep("", nrow(df))
  for(f in 1:filter_length) {
    iEmpty <- nchar(filter_states) < 1
    if(filter_ops[f] == "<") {
      iFail <- df[, filter_names[f]] < filter_values[f]
      iFail[is.na(iFail)] <- F
    } else if(filter_ops[f] == ">") {
      iFail <- df[, filter_names[f]] > filter_values[f]
      iFail[is.na(iFail)] <- F
    } else {
      print(paste0("filter input #", f, " not formatted correctly"))
    }
    filter_states[iFail & iEmpty] <- filter_names[f]
    filter_states[iFail & !iEmpty] <- paste0(filter_states[iFail & !iEmpty], ", ", filter_names[f])
  }
  return(filter_states)
}

random_round <- function(x, seed = 123, tol = 0.1, .digits = 0) { 
  set.seed(seed) 
  round(jitter(x, amount = tol), digits = .digits)
}

RMxBY_vcfParse <- function(vcf_path.) {
  # vcf_path. <- RMxBY_comp_file
  vcf_in <- read.table(vcf_path.)
  VCFhdr <- c("CHROM", "POS", "ID", "BY", "RM","QUAL", "FILTER", "INFO", "null", "FORMAT")
  colnames(vcf_in) <- VCFhdr
  vcf_in <- vcf_in[, c("CHROM", "POS", "BY", "RM")]
  # There are a few sites that have calls for two alt alleles from the RMxBY alignment. Filter out those sites.
  vcf_in <- vcf_in[grep(",", vcf_in$RM, invert=T),]
  # Some data handling requires a unique index for each position, while positions are reused among CHROMs
  chrom_indcs <- cumsum(chrom_lengths_BY)
  vcf_in$CHROM <- factor(vcf_in$CHROM)
  vcf_in$POSi <- vcf_in$POS
  # Loop to add the cumulative chromosome lengths to each CHROM. CHROM #1 requires no modificaiton
  for (chrom in 2:16) {
    #chrom<-2
    #rm(chrom)
    chromNm <- levels(vcf_in$CHROM)[chrom]
    idx <- vcf_in$CHROM==chromNm
    vcf_in[idx,"POSi"] <- vcf_in$POS[idx] + chrom_indcs[chrom-1]
  }
  return(vcf_in)
}


GenotypeFromGQ <- function(all_alleleMerge, baseThrsh = 50, 
                              naThrsh = 100, diffThrsh = 30, 
                           include_DP = T, check_alleles = T) {
  require(dplyr)
  # baseThrsh   baseline threshold - at least one call GQ must exceed this value
  # diffThrsh   difference threshold - if calls do not agree, the higher GQ 
  #             must exceed the lower GQ by this value to be called
  #---------------------------------------------------------------------------#
  # For all sites, the BY-ref-call and RM-ref-call alleles must agree
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
  
  # all_alleleMerge <- indels_merge_final %>% slice(1:1000)
  
  # Ensure all threshold variables are numeric
  baseThrsh <- as.numeric(baseThrsh)
  naThrsh <- as.numeric(naThrsh) 
  diffThrsh <- as.numeric(diffThrsh)
  
  # Include the no-call symbol "./." in the genotype levels in both callsets
  if(any(levels(all_alleleMerge$GT_BYcall) == "./.")) {
    GTlvls <- levels(all_alleleMerge$GT_BYcall)
  } else {
    GTlvls <- c(levels(all_alleleMerge$GT_BYcall), "./.")
  }
  
  # We require the genotype factors to be the same between callsets
  all_alleleMerge$GT_BYcall <- factor(all_alleleMerge$GT_BYcall, levels = GTlvls)
  all_alleleMerge$GT_RMcall <- factor(all_alleleMerge$GT_RMcall, levels = GTlvls)
  
  # Index NA rows
  iNA_BY <- is.na(all_alleleMerge$GQ_BYcall) | is.na(all_alleleMerge$GT_BYcall)
  iNA_RM <- is.na(all_alleleMerge$GQ_RMcall) | is.na(all_alleleMerge$GT_RMcall)
  
  # Convert NA to no-call
  all_alleleMerge$GT_BYcall[iNA_BY] <- "./."
  all_alleleMerge$GT_RMcall[iNA_RM] <- "./."
  
  # Index NAs for GQ and assign score of 0
  all_alleleMerge$GQ_BYcall[iNA_BY] <- 0
  all_alleleMerge$GQ_RMcall[iNA_RM] <- 0
  
  
  # All sites are initially assigned missing data values and filled in with 
  # data that passes the thresholds set
  all_alleleMerge$finalGT <- "./."
  all_alleleMerge$finalGQ <- 0
  if(include_DP) {
    all_alleleMerge$Ref_DP_final <- 0
    all_alleleMerge$Alt_DP_final <- 0
  }
  # Index same genotypes between callsets
  iNC <- all_alleleMerge$GT_BYcall == "./." & all_alleleMerge$GT_RMcall == "./."
  iSameGT <- all_alleleMerge$GT_BYcall == all_alleleMerge$GT_RMcall & !iNC
  
  # Index sites with the same call and for which the BY-ref-call has a higher GQ value 
  iBYoverRM <- iSameGT & all_alleleMerge$GQ_BYcall >= all_alleleMerge$GQ_RMcall
  
  # Likewise with the RM-ref-call set
  iRMoverBY <- iSameGT & all_alleleMerge$GQ_RMcall > all_alleleMerge$GQ_BYcall
  
  # Assign values by indicies from above
  all_alleleMerge$finalGT[iBYoverRM] <- as.character(all_alleleMerge$GT_BYcall[iBYoverRM])
  all_alleleMerge$finalGT[iRMoverBY] <- as.character(all_alleleMerge$GT_RMcall[iRMoverBY])
  all_alleleMerge$finalGQ[iBYoverRM] <- all_alleleMerge$GQ_BYcall[iBYoverRM]
  all_alleleMerge$finalGQ[iRMoverBY] <- all_alleleMerge$GQ_RMcall[iRMoverBY]
  
  if(include_DP) {
    all_alleleMerge$Ref_DP_final[iSameGT] <- random_round(rowSums(all_alleleMerge[iSameGT, c("Ref_DP_BYcall", "Ref_DP_RMcall")])/2)
    all_alleleMerge$Alt_DP_final[iSameGT] <- random_round(rowSums(all_alleleMerge[iSameGT, c("Alt_DP_BYcall", "Alt_DP_RMcall")])/2)
  }
  # Assign genotypes for discordant sites. The genotype for each site is assigned
  # as either the BY-ref- or RM-ref-derived call, depending on which has a higher
  # GQ value.
  # Index sites for which the BY call will be used with the logic:
  # genotypes do not agree -and- 
  # # the BY call is not NA -and- 
  # # # either RM call is not NA and BY GQ exceeds RM by the diffThrsh -or-
  # # # RM is NA and BY GQ exceeds the naThrsh
  iBYcall <- !iSameGT & 
    !iNA_BY & 
    ((!iNA_RM & 
        all_alleleMerge$GQ_BYcall > all_alleleMerge$GQ_RMcall + diffThrsh) | 
       (iNA_RM & 
          all_alleleMerge$GQ_BYcall >= naThrsh))
  
  iRMcall <- !iSameGT & 
    !iNA_RM &
    ((!iNA_BY & 
        all_alleleMerge$GQ_RMcall > all_alleleMerge$GQ_BYcall + diffThrsh) |
       (iNA_BY & 
          all_alleleMerge$GQ_RMcall >= naThrsh))
  
  all_alleleMerge$finalGT[iBYcall] <- as.character(all_alleleMerge$GT_BYcall[iBYcall])
  all_alleleMerge$finalGQ[iBYcall] <- all_alleleMerge$GQ_BYcall[iBYcall]
  all_alleleMerge$finalGT[iRMcall] <- as.character(all_alleleMerge$GT_RMcall[iRMcall])
  all_alleleMerge$finalGQ[iRMcall] <- all_alleleMerge$GQ_RMcall[iRMcall]
  
  if(include_DP) {
    all_alleleMerge$Ref_DP_final[iBYcall] <- all_alleleMerge$Ref_DP_BYcall[iBYcall]
    all_alleleMerge$Alt_DP_final[iBYcall] <- all_alleleMerge$Alt_DP_BYcall[iBYcall]
    all_alleleMerge$Ref_DP_final[iRMcall] <- all_alleleMerge$Ref_DP_RMcall[iRMcall]
    all_alleleMerge$Alt_DP_final[iRMcall] <- all_alleleMerge$Alt_DP_RMcall[iRMcall]
  }
  if(check_alleles) {
  # Revert sites for which alleles do not agree within the set that have received a final call
    all_alleleMerge$REF_BYcall[is.na(all_alleleMerge$REF_BYcall)] <- "."
    all_alleleMerge$ALT_BYcall[is.na(all_alleleMerge$ALT_BYcall)] <- "."
    all_alleleMerge$REF_RMcall[is.na(all_alleleMerge$REF_RMcall)] <- "."
    all_alleleMerge$ALT_RMcall[is.na(all_alleleMerge$ALT_RMcall)] <- "."
    
    iDiffAllele <-  all_alleleMerge$REF_BYcall != "." & 
      all_alleleMerge$REF_RMcall != "." &
      !(all_alleleMerge$REF_BYcall == all_alleleMerge$REF_RMcall &
          all_alleleMerge$ALT_BYcall == all_alleleMerge$ALT_RMcall)
    # all_alleleMerge[iDiffAllele, ]
    
    all_alleleMerge$finalGT[iDiffAllele] <- "./."
    all_alleleMerge$finalGQ[iDiffAllele] <- 0
    if(include_DP) {
      all_alleleMerge$Ref_DP_final[iDiffAllele] <- 0
      all_alleleMerge$Alt_DP_final[iDiffAllele] <- 0
    }
  }
  # Revert sites not exceeding baseThrsh
  iBYunderBase <- all_alleleMerge$GQ_BYcall < baseThrsh
  iBYunderBase[is.na(iBYunderBase)] <- T
  iRMunderBase <- all_alleleMerge$GQ_RMcall < baseThrsh
  iRMunderBase[is.na(iRMunderBase)] <- T
  iUnderBase <- iBYunderBase & iRMunderBase
  
  all_alleleMerge$finalGT[iUnderBase] <- "./."
  all_alleleMerge$finalGQ[iUnderBase] <- 0
  if(include_DP) {
    all_alleleMerge$Ref_DP_final[iUnderBase] <- 0
    all_alleleMerge$Alt_DP_final[iUnderBase] <- 0
  }
  all_alleleMerge$finalGT <- factor(all_alleleMerge$finalGT, levels = GTlvls)
  return(all_alleleMerge)
}

MarkDubiousHets <- function(SNP_df, .alpha = 0.05, known_hets = "empirical", bin_size = 0.1) {
  # bin_size is in log10 scale. A value of 0.1 yields approximately 30 bins
  # SNP_df <- indels_merge_finalGT
  # GTstats <- all_wide_GTstats
  if(length(known_hets) > 2 & is.numeric(known_hets)) {
    # We can also provide known heterozygous sites to the function instead.
    anc_het_POSi <- known_hets
    known_het_sites <- SNP_df %>% select(POSi, ID, Rep, GT, Sum_DP_final, f_Alt) %>%
      filter(POSi %in% known_hets, GT == "0/1", Rep == "00") %>% arrange(ID, POSi)
  } else if(known_hets == "empirical") {
    # We can obtain known heterozygous sites from the founder data. If at least
    # 7 founders are called and all calls are heterozygous, we can be reasonably
    # sure that they are true ancestral heterozygous sites.
    anc_GTstats <- SNP_df %>% site_genotype_stats(group = "anc")
    anc_het_POSi <- anc_GTstats %>% 
      filter(fHet == 1, nHet >= 7) %>% 
      pull(POSi)
    known_het_sites <- SNP_df %>% select(POSi, ID, Rep, GT, Sum_DP_final, f_Alt) %>%
      filter(POSi %in% anc_het_POSi, GT == "0/1", Rep == "00") %>% arrange(ID, POSi)
  } else {
    print("No known heterozygous sites available")
  }

  bin_by_DP <- function(sites, DP_col = "Sum_DP_final", return_bin = F, bin_input = "from_data") {
    # Get range in read depth at Het sites, then form bins
    # sites <- SNP_df
    if(is.character(bin_input)) {
      if(bin_input ==  "from_data") {
        temp_range <- range(log10(sites$Sum_DP_final))
        bin_range = c(floor(log10(6)*100)/100, ceiling(temp_range[2]*100)/100)
        log_bins = seq(from = bin_range[1], to = bin_range[2], by = bin_size)
        log_bins <- c(0, log_bins[-c(1,length(log_bins))])
        bin_max <- round(max(log10(sites$Sum_DP_final)) + 3, 2)
        log_bins_df <- data.frame(bin = log_bins, upper_bin = c(log_bins[-1], bin_max))
        log_bins_df$label <- paste0(log_bins_df$bin, " - ", log_bins_df$upper_bin)
      } else {
        print("bin input not recognized")
      }
    } else {
      log_bins_df <- bin_input
    }
    if(return_bin) {
      if(bin_input == "from_data") {
        return(log_bins_df)
      } else {
        print("No data for bin generation")
      }
    } else {
    sites$log_DP <- log10(sites[, DP_col])
    # For each bin range, label sites with bin membership
    sites$bin <- 0
    for(b in 1:nrow(log_bins_df)) {
      # b = 2
      i_bin <- sites$log_DP >= log_bins_df$bin[b] & 
        sites$log_DP < log_bins_df$upper_bin[b]
      sites$bin[i_bin] <- log_bins_df$bin[b]
    }
    
    # Bins with fewer than 100 sites should be merged into neighboring groups
    small_bins <- sites %>% 
      count(bin) %>% filter(n < 100)
    bin_counts <- sites %>% 
      count(bin) 
    if(nrow(small_bins) != 0) {
      med_bin <- median(sites$bin)
      small_bins_lo <- small_bins %>% filter(bin < med_bin)
      if(nrow(small_bins_lo) > 1) {
        sbl_cont <- sum(diff(which(log_bins_df$bin %in% small_bins_lo$bin)) != 1) == 0
      } else if(nrow(small_bins_lo) == 1) {
        sbl_cont <- T
      } else {
        sbl_cont <- F
        print("Low small bins not present")
      }
      
      small_bins_hi <- small_bins %>% filter(bin > med_bin)
      small_bins_hi <- small_bins_hi %>% arrange(desc(bin))
      if(nrow(small_bins_hi) > 1) {
        sbh_cont <- sum(diff(which(log_bins_df$bin %in% small_bins_hi$bin)) != 1) == 0
      } else if(nrow(small_bins_hi) == 1) {
        sbh_cont <- T
      } else {
        sbh_cont <- F
        print("High small bins not present")
      }
      
      if(nrow(small_bins_lo) > 0 & sbl_cont) {
        max_bin_lo <- max(small_bins_lo$bin)
        to_bin <- min(sites$bin[sites$bin > max_bin_lo])
        i_merge <- sites$bin <= max_bin_lo
        sites$bin[i_merge] <- to_bin
      } else if(nrow(small_bins_lo) > 0 & !sbl_cont) {
        for(i in seq_along(small_bins_lo$bin)) {
          bin_lo <- small_bins_lo$bin[i]
          to_bin <- min(sites$bin[sites$bin > bin_lo])
          i_merge <- sites$bin == bin_lo
          sites$bin[i_merge] <- to_bin
        }
      } else {
        print("Low small bin not processed")
      }
      
      if(nrow(small_bins_hi) > 0 & sbh_cont) {
        min_bin_hi <- min(small_bins_hi$bin)
        to_bin <- max(sites$bin[sites$bin < min_bin_hi])
        i_merge <- sites$bin >= min_bin_hi
        sites$bin[i_merge] <- to_bin
      } else if(nrow(small_bins_hi) > 0 & !sbh_cont) {
        for(i in seq_along(small_bins_hi$bin)) {
          bin_hi <- small_bins_hi$bin[i]
          to_bin <- min(sites$bin[sites$bin < bin_lo])
          i_merge <- sites$bin == bin_hi
          sites$bin[i_merge] <- to_bin
        }
      } else {
        print("High small bin not processed")
      }
    }
    # Include bin labels for plotting
    final_bins <- unique(sites$bin) %>% sort()
    bin_names <- paste0(final_bins, " - ", c(final_bins[-1], "inf"))
    names(bin_names) <- final_bins
    sites$bin_name <- factor(sites$bin, labels = bin_names)
    return(sites)
    }
  }
  
  known_het_sites <- known_het_sites %>% bin_by_DP()
  known_het_bins <- known_het_sites %>% bin_by_DP(return_bin = T)
  final_bins <- known_het_bins$bin
  bin_names <- known_het_bins$label
  # Calculate mean and variance of the alt read frequency for each bin

  SNP_df <- SNP_df %>% bin_by_DP(bin_input = known_het_bins)
  
  mean_alt_freq_known_het <- known_het_sites %>% group_by(bin, bin_name) %>% 
    summarize(f_alt = mean(f_Alt), var_alt = var(f_Alt))
  
  # Use variance in alt allele frequency for each bin to calculate
  # standard deviation of the null normal distribution
  # Expected alt allele frequency is set to 0.5
  .sigma <- mean_alt_freq_known_het$var_alt %>% sqrt()
  mean_alt_freq_known_het$low_cut <- qnorm(p = .alpha/2, mean = 0.5, sd = .sigma)
  mean_alt_freq_known_het$high_cut <- qnorm(p = 1 - .alpha/2, mean = 0.5, sd = .sigma)
  
  # Form logical vector of sites outside null distribution threshold .alpha
  cut_out <- rep(F, nrow(SNP_df))
  GT_split <- colsplit(SNP_df$GT, "/", names = c("A_1", "A_2"))
  for(b in 1:nrow(mean_alt_freq_known_het)) {
    # b = 1
    i_bin <- SNP_df$bin == mean_alt_freq_known_het$bin[b]
    cuts <- c(mean_alt_freq_known_het$low_cut[b], mean_alt_freq_known_het$high_cut[b])
    # Index current bin, Het sits, outside .alpha, and not in the training set
    i_cut <- i_bin & GT_split$A_1 != GT_split$A_2 & # in bin and het
      (SNP_df$f_Alt <= cuts[1] | SNP_df$f_Alt >= cuts[2]) & # outside of allele fraction CIs
      !(SNP_df$POSi %in% anc_het_POSi) # not in training set
    cut_out[i_cut] <- T
  }
  i_zeroAlt_het <- SNP_df$f_Alt == 0 & SNP_df$GT == "0/1" & !(SNP_df$POSi %in% anc_het_POSi)
  cut_out[i_zeroAlt_het] <- T
  return(cut_out)
}

FormatSNPsMerge <- function(SNPs_finalGT_df) {
  if(sum(colnames(SNPs_finalGT_df) == "Sum_DP_final") == 0) {
  SNPs_finalGT_df$Sum_DP_final <- SNPs_finalGT_df$Ref_DP_final + SNPs_finalGT_df$Alt_DP_final
  }
  SNPs_finalGT_df$Line <- droplevels(SNPs_finalGT_df$Line)
  SNPs_finalGT_df$ID <- droplevels(SNPs_finalGT_df$ID)
  SNPs_finalGT_df$Tx_name <- Recode_Tx(SNPs_finalGT_df$Tx)
  SNPs_finalGT_df$Tx_ID <- Recode_Tx_ID(SNPs_finalGT_df$Tx, tx_type = "Tx")
  GTGQ_cols <- grep("finalG", colnames(SNPs_finalGT_df))
  colnames(SNPs_finalGT_df)[GTGQ_cols] <- c("GT", "GQ")
  if(sum(colnames(SNPs_finalGT_df) == "existing_SNP") == 0) {
    SNPs_finalGT_df$existing_SNP <- F
    i_eSNP <- SNPs_finalGT_df$POSi %in% RMxBY_comp_SNPs_POSi
    SNPs_finalGT_df$existing_SNP[i_eSNP] <- T
  }
  if(sum(colnames(SNPs_finalGT_df) == "f_Alt") == 0) {
    SNPs_finalGT_df$f_Alt <- SNPs_finalGT_df$Alt_DP_final/SNPs_finalGT_df$Sum_DP_final
  }
  return(SNPs_finalGT_df)
}

InferMissingFounders <- function(SNPs_finalGT_df, missing_founders = noAncestor, het_support = 4, hom_support = 7) {
  # SNPs_finalGT_df <- SNPs_merge_finalGT
  missing_founders <- as.character(missing_founders)
  Anc_add <- data.frame(NULL)
  for(i in seq_along(missing_founders)) {
    # i = 2
    # Get clones from line missing founder
    SNPs_noAnc <- SNPs_finalGT_df %>% 
      filter(Line == missing_founders[i], Rep != "00", !Cut)
    noAnc_line <- missing_founders[i]
    # Get counts of each genotype across these clones
    evo_GTstats <- SNPs_noAnc %>% site_genotype_stats(group = "evo")
    # Index columns with position info
    pos_cols <- c(grep("CHROM", colnames(evo_GTstats)), grep("POS", colnames(evo_GTstats)))
    # Get genotype counts across all existing founder clones
    anc_GTstats <- SNPs_finalGT_df %>%
      filter(!Cut) %>%
      site_genotype_stats(group = "anc")
    all_GTstats <- merge(anc_GTstats, evo_GTstats, all = T,
                           by = c("CHROM", "POS", "POSi"), 
                           suffixes = c("_anc", "_evo")) %>% 
      filter(!is.na(fHet_anc), !is.na(fHet_evo))
    
    # all_GTstats <- cbind(evo_GTstats[, pos_cols], anc_GTstats, 
    #                      evo_GTstats %>% select(contains("evo")))
    
    # Get consensus genotypes for each position and index, construct list of indicies
    i_het <- all_GTstats$fHet_anc == 1 & 
      all_GTstats$nHet_evo > het_support & all_GTstats$nNA_evo < 10
    i_Ref <- all_GTstats$fRef_anc == 1 & 
      all_GTstats$nRef_evo > hom_support & all_GTstats$nNA_evo < 10
    i_Alt <- all_GTstats$fAlt_anc == 1 & 
      all_GTstats$nAlt_evo > hom_support & all_GTstats$nNA_evo < 10
    gt_list <- list(i_het, i_Ref, i_Alt)
    names(gt_list) <- c("0/1", "0/0", "1/1")
    
    # all_GTstats$anc_GT <-  "./."
    # Loop over list of positions with each genotype. For each genotype, pull 
    # end-point clone rows with matching genotypes and positions. From this 
    # set of clones and positions, randomly choose rows without replacement, 
    # remove duplicate positions, and add to founder. Remove rows that have a
    # match in the founder. Repeat until there are not any unmatched positions 
    # in the founder. Repeat for all genotypes.
    anc_rows <- data.frame(NULL)
    for(gt in names(gt_list)) {
      # gt <- names(gt_list)[1]
      anc_gt_rows <- data.frame(NULL)
      # all_GTstats$anc_GT[gt_list[[gt]]] <- gt
      gt_POSi <- all_GTstats %>% 
        filter(gt_list[[gt]]) %>% pull(POSi)
      l <- length(gt_POSi)
      anc_sample <- SNPs_noAnc %>% 
        filter(GT == gt & POSi %in% gt_POSi)
      while(l > 0) {
        gt_random <- anc_sample %>% 
          sample_n(l, replace = F) %>% 
          distinct(POSi, .keep_all = T)
        anc_gt_rows <- rbind(anc_gt_rows, gt_random)
        anc_sample <- anc_sample %>% 
          filter(!POSi %in% anc_gt_rows$POSi)
        n_done <- nrow(gt_random)
        l <- l - n_done
      }
      anc_rows <- rbind(anc_rows, anc_gt_rows)
    }
    anc_rows <- anc_rows %>% distinct(POSi, .keep_all = T)
    anc_rows$ID <- paste0(noAnc_line, "00")
    anc_rows$Rep<- "00"
    Anc_add <- rbind(Anc_add, anc_rows)
  }
  
  # Anc_add$Rep<- "00"
  # Anc_add$ID <- paste0(Anc_add$Line, "00")
  
  noAnc_merge_finalGT <- rbind(Anc_add %>% 
                                 filter(!(GT != "0/1" & existing_SNP)), 
                               SNPs_finalGT_df %>% 
                                 filter(Line %in% missing_founders, Rep != "00")) %>%
                            arrange(ID, POSi)
  
  # noAnc_merge_finalGT %>% site_genotype_stats()
  
  # noAnc_merge_finalGT$Tx <- factor(noAnc_merge_finalGT$Tx)
  # noAnc_merge_finalGT$Line <- factor(noAnc_merge_finalGT$Line)
  # noAnc_merge_finalGT$Rep <- factor(noAnc_merge_finalGT$Rep)
  # noAnc_merge_finalGT$ID <- factor(noAnc_merge_finalGT$ID)
  # noAnc_merge_finalGT$GT <- factor(noAnc_merge_finalGT$GT)
  
  SNPs_finalGT_df <- rbind(noAnc_merge_finalGT, 
                           SNPs_finalGT_df %>% filter(!Line %in% missing_founders)) 
  
  SNPs_finalGT_df$ID <- as.character(SNPs_finalGT_df$ID)
  SNPs_finalGT_df$Line <- as.character(SNPs_finalGT_df$Line)
  SNPs_finalGT_df$Rep <- as.character(SNPs_finalGT_df$Rep)

  SNPs_finalGT_df$Tx <- factor(SNPs_finalGT_df$Tx)
  SNPs_finalGT_df$Line <- factor(SNPs_finalGT_df$Line)
  SNPs_finalGT_df$Rep <- factor(SNPs_finalGT_df$Rep)
  SNPs_finalGT_df$ID <- factor(SNPs_finalGT_df$ID)
  SNPs_finalGT_df$GT <- factor(SNPs_finalGT_df$GT)
  
  SNPs_finalGT_df <- SNPs_finalGT_df %>%
    arrange(ID, POSi)
  
  return(SNPs_finalGT_df)
}

anc_GT_fltr <- function(all_vcf_in, anc_GT = "0/1", min_anc_GQ = 0, min_evo_GQ = 0) {
  ## Takes in a dataframe containing genotypes for founder and end-point clones 
  ## and filters end-point sites to only those that were of a set genotype (anc_GT) in the founder. 
  ## Optional: First filters foudner clones to quality threshold (GQ >= "ancQual") and then
  ## evolved clones (GenQual >= "evoQual") at sites where the founder is of the specified genotype
  require(dplyr)
  # all_vcf_in <- noAnc_merge_finalGT
  all_vcf_in$Line <- factor(all_vcf_in$Line)
  all_vcf_in$ID <- factor(all_vcf_in$ID)
  all_vcf_in <- all_vcf_in %>%
    arrange(ID, POSi)
  for (l in 1:length(levels(all_vcf_in$Line))) {
    # l=1
    line_GTfltr <- all_vcf_in %>% 
      filter(Line == levels(all_vcf_in$Line)[l]) 
    
    if(sum(line_GTfltr$Rep == "00") == 0) {
      anc_sites <- line_GTfltr %>% 
        filter(GT == anc_GT) %>% 
        distinct(POSi) 
    } else {
      anc_sites <- line_GTfltr %>% 
        filter(Rep == "00" & GT == anc_GT)
      if (min_anc_GQ == 0) {
        anc_sites <- anc_sites %>%
          select(POSi)
      } else {
        anc_sites <- anc_sites %>%
          filter(GQ >= min_anc_GQ) %>%
          select(POSi)
      }
    }
    if (min_evo_GQ == 0) {
      evo_fltrd <- line_GTfltr %>%
        filter(POSi %in% anc_sites$POSi)
    } else {
      evo_fltrd <- line_GTfltr %>%
        filter(POSi %in% anc_sites$POSi) %>%
        filter(GQ >= min_evo_GQ)
    }

    if (l == 1) {
      evo_out <- evo_fltrd
    } else {
      evo_out <- rbind(evo_out, evo_fltrd)
    }
    print(levels(all_vcf_in$Line)[l])
  }
  # evo_out$chrom_n <- as.numeric(evo_out$CHROM)
  # evo_out$CHROM <- factor(str_pad(evo_out$chrom_n, width = 2, pad = "0"))
  evo_out$GT <- droplevels(evo_out$GT)
  evo_out$ID <- droplevels(evo_out$ID)
  evo_out$isHom <- as.numeric(evo_out$GT != "0/1")
  return(evo_out)
}

is_outlier <- function(x) {
  ## Returns a logical vector with true values for data points that are more extreme than 
  ## 1.5 * IQR outside the 0.25 and 0.75 quantiles
  return(x < quantile(x, 0.25) - 1.5 * IQR(x) | x > quantile(x, 0.75) + 1.5 * IQR(x))
}

errorFromPhylo <- function(vcf_df, flsHom_support = 4, flsHet_support = 6, output_POSi = F) {
  # vcf_df <- SNPs_merge_finalGT 
  vcf_df$Line <- factor(vcf_df$Line)
  if(sum(levels(vcf_df$GT) == "./.") == 0) {
    vcf_df$GT <- factor(vcf_df$GT, levels = c(levels(vcf_df$GT), "./."))
  }
  POSi_out <- data.frame(NULL)
  df_out <- data.frame(NULL)
  for(fg in levels(vcf_df$Line)) {
    # fg <- "N_A"
    groupVcf <- vcf_df %>% filter(Line == fg)
    groupVcf <- groupVcf %>% 
                  group_by(Rep) %>% 
                  filter(!duplicated(POSi))
    
    
    evo_GT_wide <- groupVcf %>% 
                      # filter(Rep != "00" & GT == "het") %>% 
                      filter(!is.na(POSi)) %>%
                      ungroup() %>%
                      select(CHROM, POS, POSi, ID, GT) %>% 
                      pivot_wider(names_from = ID, values_from = GT)
    colnames(evo_GT_wide)[4] <- "anc_GT"
    evo_GT_wide$anc_GT[is.na(evo_GT_wide$anc_GT)] <- "./."
    nAnc_sites <- sum(evo_GT_wide$anc_GT != "./.")
    n_evo <- ncol(evo_GT_wide) - 4
    n_het <- apply(evo_GT_wide[, -c(1:4)], 1, 
                               function(x) sum(x == "0/1", na.rm = T))
    n_Ref <- apply(evo_GT_wide[, -c(1:4)], 1, 
                               function(x) sum(x == "0/0", na.rm = T))
    n_Alt <- apply(evo_GT_wide[, -c(1:4)], 1, 
                               function(x) sum(x == "1/1", na.rm = T))
    n_na <- apply(evo_GT_wide[, -c(1:4)], 1, 
                               function(x) sum(is.na(x)) + sum(x == "./.", na.rm = T))
    evo_GT_wide <- cbind(evo_GT_wide, nRef = n_Ref, 
                         nHet = n_het, nAlt = n_Alt, nNA = n_na)
    n_Het_q <- sum(evo_GT_wide$anc_GT != "./." & evo_GT_wide$nHet >= flsHom_support)
    iRef <- evo_GT_wide$anc_GT == "0/0" & evo_GT_wide$nHet >= flsHom_support
    n_Ref_hom <- sum(iRef, na.rm = T)
    iHet <- evo_GT_wide$anc_GT == "0/1" & evo_GT_wide$nHet >= flsHom_support
    n_het <- sum(iHet, na.rm = T)
    iAlt <- evo_GT_wide$anc_GT == "1/1" & evo_GT_wide$nHet >= flsHom_support
    n_Alt_hom <- sum(iAlt, na.rm = T)
    f_Ref <- n_Ref_hom/(n_Ref_hom + n_Alt_hom)
    
    i_Hom_q <- (evo_GT_wide$nRef >= flsHet_support | evo_GT_wide$nAlt >= flsHet_support) & 
      evo_GT_wide$nHet == 0
    n_Hom_q <- sum(evo_GT_wide$anc_GT != "./." & i_Hom_q)
    iRef_F_het <- evo_GT_wide$anc_GT == "0/1" & evo_GT_wide$nRef >= flsHet_support & evo_GT_wide$nHet == 0
    n_Ref_F_het <- sum(iRef_F_het, na.rm = T)
    iAlt_F_het <- evo_GT_wide$anc_GT == "0/1" & evo_GT_wide$nAlt >= flsHet_support & evo_GT_wide$nHet == 0
    n_Alt_F_het <- sum(iAlt_F_het, na.rm = T)
    # i_sample_FN <- (evo_GT_wide$nRef >= flsHet_support | evo_GT_wide$nAlt >= flsHet_support) & evo_GT_wide$nHet == 0
    # n_sample_FN <- ifelse(sum(i_sample_FN, na.rm = T) == 0, 1, sum(i_sample_FN, na.rm = T))
    if(output_POSi == T) {
      POSi_df <- cbind(Line = fg, 
                       evo_GT_wide[, c("CHROM", "POS", "POSi", "anc_GT",
                                       "nRef", "nHet", "nAlt", "nNA")], 
                       is_error = (iRef | iAlt))
      POSi_out <- rbind(POSi_out, POSi_df)
      print(fg)
    } else {
      fg_df <- data.frame(Line = fg, nHet = n_het, 
                          n_Het_q = n_Het_q,
                          nRef = n_Ref_hom, Ref_rate = n_Ref_hom/n_Het_q, 
                          nAlt = n_Alt_hom, Alt_rate = n_Alt_hom/n_Het_q, 
                          n_F_Hom = n_Ref_hom + n_Alt_hom, 
                          F_Hom_rate = (n_Ref_hom + n_Alt_hom)/n_Het_q, 
                          Ref_bias = f_Ref, 
                          n_Hom_q = n_Hom_q,
                          nRef_F_Het = n_Ref_F_het, Ref_F_Het_rate = n_Ref_F_het/(n_Hom_q),
                          nAlt_F_Het = n_Alt_F_het, Alt_F_Het_rate = n_Alt_F_het/(n_Hom_q),
                          n_F_Het = n_Ref_F_het + n_Alt_F_het, 
                          F_Het_rate = (n_Ref_F_het + n_Alt_F_het)/(n_Hom_q))
      df_out <- rbind(df_out, fg_df)
      print(fg)
    }
  }
  if(output_POSi == T) {
    POSi_out$Line <- factor(POSi_out$Line)
    return(POSi_out)
  } else {
    df_out$Line <- factor(df_out$Line)
    return(df_out)
  }
}

EstDataBounds <- function(anc_het_df, rm_noData = T) {
  require(pammtools)
  # Takes the "all_ancHet" dataframe and estimates the boundaries
  # of LOH regions and no data regions
  # Output is a table of start and end positions of contiguous regions
  # with t_i = -1 for no data and t_i = as.numeric(GT) for called genotypes
  
  # anc_het_df <- LOH_SNPs %>% filter(Rep != "00", Line == "F_C") 
  
  anc_het_df$ID <- factor(anc_het_df$ID)
  anc_het_df$GT <- factor(anc_het_df$GT)
  # Vector of positions called across dataset
  all_POSi <- anc_het_df %>% 
    select(POSi, CHROM, POS) %>% 
    distinct() %>% 
    arrange(POSi)
  t_pos_df <- data.frame(NULL)
  for (id in levels(anc_het_df$ID)) {
    # Subset by clone
    # id = "F_C07"
    id_ancHet <- anc_het_df %>% filter(ID == id) %>% 
      select(ID, POSi, CHROM, POS, GT) %>% 
      arrange(POSi)
    # Merge clone sites with positions vector
    id_SNPs <- merge(all_POSi, id_ancHet, 
                      by = c("POSi", "CHROM", "POS"), 
                      sort = T, all = T) 
    # Add columns classifying sites as: 
    # t_i = -1    no data
    # t_i != -1   GT levels
    id_SNPs <- id_SNPs %>% 
      mutate(a_i = as.numeric(GT), 
             d_i = as.numeric(!is.na(GT)))
    id_SNPs$ID <- id
    id_SNPs$t_i <- -1
    id_SNPs$t_i[id_SNPs$d_i == 1] <- as.numeric(id_SNPs$GT[id_SNPs$d_i == 1])
    id_SNPs$CHROM <- factor(id_SNPs$CHROM)
    id_SNPs$CHROM <- droplevels(id_SNPs$CHROM)
    # Detect runs of t_i and create a dataframe of positions
    if (rm_noData == T) {
      id_SNPs <- id_SNPs %>% filter(t_i != -1)
    }
    t_runs_df <- data.frame(NULL)
    for (ch in levels(id_SNPs$CHROM)) {
      # ch = "09"
      t_runs <- rle(id_SNPs$t_i[id_SNPs$CHROM == ch])
      t_runs_ch <- data.frame(value = t_runs$values, 
                              length = t_runs$lengths,
                              CHROM = ch)
      t_runs_df <- rbind(t_runs_df, t_runs_ch)
    }
    t_runs_df$CHROM <- factor(str_pad(t_runs_df$CHROM, width = 2, pad = "0"))
    t_runs_df$start <- cumsum(t_runs_df$length) - t_runs_df$length + 1
    t_runs_df$end <-  cumsum(t_runs_df$length)
    t_runs_df$start_POSi <- id_SNPs$POSi[t_runs_df$start]
    t_runs_df$end_POSi <- id_SNPs$POSi[t_runs_df$end]
    # Estimate boundaries of type events
    t_runs_df$est_start <- 0
    t_runs_df$est_end <- 0
    for (i in 1:nrow(t_runs_df)) {
      # i = 24
      obsv_start <- t_runs_df$start_POSi[i]
      obsv_end <- t_runs_df$end_POSi[i]
      if (i == 1) {
        t_runs_df$est_start[i] <- 1
      } else {
        if(t_runs_df$CHROM[i] != t_runs_df$CHROM[i - 1]) {
          t_runs_df$est_start[i] <- chrom_bound_BY %>% filter(CHROM == t_runs_df$CHROM[i]) %>% pull(Start)
        } else {
          t_runs_df$est_start[i] <- ceiling((obsv_start + t_runs_df$end_POSi[i - 1]) / 2)
        }
      } 
      if (i == nrow(t_runs_df)) {
        t_runs_df$est_end[i] <- chrom_bound_BY$End[16]
      } else {
        if(t_runs_df$CHROM[i] != t_runs_df$CHROM[i + 1]) {
          t_runs_df$est_end[i] <- chrom_bound_BY %>% filter(CHROM == t_runs_df$CHROM[i]) %>% pull(End)
        } else {
        t_runs_df$est_end[i] <- floor((obsv_end + t_runs_df$start_POSi[i + 1] - 1) / 2)
        }
      }
    }
    t_runs_df$est_length <- t_runs_df$est_end - t_runs_df$est_start + 1
    t_runs_df <- t_runs_df %>% filter(value != -1) 
    t_runs_df$GT <- levels(anc_het_df$GT)[t_runs_df$value]
    t_runs_df <- t_runs_df %>% mutate(ID = id)
    t_pos_df <- rbind(t_pos_df, t_runs_df)
    # print(id)
  } 
  t_pos_df$CHROM <- factor(t_pos_df$CHROM)
  t_pos_df$ID <- factor(t_pos_df$ID)
  t_pos_df$Tx <- factor(substr(t_pos_df$ID, 1, 1))
  t_pos_df$Tx_name <- Recode(t_pos_df$Tx, "WT" = "N", "Cas9" = "H", "Drive" = "F")
  t_pos_df$Line <- factor(substr(t_pos_df$ID, 1, 3))
  t_pos_df$Rep <- factor(substr(t_pos_df$ID, 4, 5))
  return(t_pos_df)
}

MarkLOHcomplexes <- function(dataBounds_df, gap = 10000) {
  # dataBounds_df <- all_LOHbounds
  dataBounds_df$ID <- factor(dataBounds_df$ID)
  
  out_df <- data.frame(NULL)
  for (id in levels(dataBounds_df$ID)) {
    # id = "N_A05"
    id_df <- subset(dataBounds_df, ID == id)
    id_df$LOH_cmplx <- 0
    id_df$LOH_k <- 0
    for (ch in levels(id_df$CHROM)) {
      # ch = "09"
      chrom_df <- id_df %>% 
        filter(CHROM == ch) %>%
        filter(GT != "0/1" & GT != "./.")
      if (nrow(chrom_df) >= 2) {
        for(l in 1:(nrow(chrom_df) - 1)) {
          # l <- 1
          l_start <- chrom_df$est_start[l]
          q_start <- chrom_df$est_end[l]
          q_end <- chrom_df$est_start[l + 1]
          l_end <- chrom_df$est_end[l + 1]
          if (q_end - q_start <= gap) {
            id_df$LOH_cmplx[id_df$est_start >= l_start & id_df$est_end <= l_end] <- 1
          }
        }
      }
    }
    k <- 1
    for (r in 1:(nrow(id_df) - 1)) {
      # r <- 27
      if (id_df$LOH_cmplx[r] == 1) {
        id_df$LOH_k[r] <- k
      } 
      if (id_df$LOH_cmplx[r] == 1 & (id_df$LOH_cmplx[r + 1] == 0 | 
           id_df$CHROM[r] != id_df$CHROM[r + 1])) {
        k <- k + 1
      }
    } 
    if(id_df$LOH_cmplx[nrow(id_df)] == 1) {
      id_df$LOH_k[nrow(id_df)] <- id_df$LOH_k[nrow(id_df) - 1]
    }
    out_df <- rbind(out_df, id_df)
  }
  return(out_df)
}

MergeComplexLOHs <- function(markedBounds_df) {
  # markedBounds_df <- all_LOHbounds %>% MarkLOHcomplexes()
  # markedBounds_df <- markedBounds_df2 %>% filter(GT != "0/1")
  end_cols <- c("end", "end_POSi", "est_end")
  for(id in levels(markedBounds_df$ID)) {
    # id <- "N_A05"
    i_id <- which(levels(markedBounds_df$ID) == id)
    # id = levels(markedBounds_df$ID)[25]
    id_max_k <- markedBounds_df %>% filter(ID == id) %>% pull(LOH_k) %>% max()
    if(id_max_k > 0) {
      for(k in 1:id_max_k) {
        # k = 2
        # index the rows to include in the complex LOH, retain first row and
        # transfer LOH end position info from last row, recalculate length
        i_complex <- which(markedBounds_df$ID == id & markedBounds_df$LOH_k == k)
        i_biggest <- i_complex[which.max(markedBounds_df$length[i_complex])]
        markedBounds_df[min(i_complex), end_cols] <- markedBounds_df[max(i_complex), end_cols]
        markedBounds_df$GT[min(i_complex)] <- markedBounds_df$GT[i_biggest]
        i_length <- markedBounds_df$end[max(i_complex)] - markedBounds_df$start[min(i_complex)]
        est_length <- markedBounds_df$est_end[max(i_complex)] - markedBounds_df$est_start[min(i_complex)]
        markedBounds_df$length[min(i_complex)] <- i_length
        markedBounds_df$est_length[min(i_complex)] <- est_length
        markedBounds_df <- markedBounds_df[-i_complex[2:length(i_complex)], ]
      }
    }
    if(i_id %% 10 == 0) {
      print(paste0(round(i_id/length(levels(markedBounds_df$ID)), 2)*100, "% done"))
    }
  }
  return(markedBounds_df)
}

MarkErrorLOH <- function(markedBounds_df, error_rate) {
  # markedBounds_df <- all_LOHbounds_merge_EC
  # error_rate <- overall_F_Hom_rate
  rownames(markedBounds_df) <- NULL
  markedBounds_df$ID <- factor(markedBounds_df$ID)
  
  id_nSites <- markedBounds_df %>% 
    group_by(ID) %>% 
    summarise(nSites = sum(length))
  
  nFP_estimate <- data.frame(ID = factor(id_nSites$ID), 
                             nFP = ceiling(id_nSites$nSites * error_rate))
  markedBounds_df$is_error <- FALSE
  singles_df <- markedBounds_df %>% filter(length == 1, GT != "0/1")
  all_others_df <- markedBounds_df %>% filter(!(length == 1 & GT != "0/1"))
  # numError <- data.frame(ID = levels(markedBounds_df$ID), nError = nFP_estimate$nFP)
  singles_mrkErr <- data.frame(NULL)
  # nErrorLOH <- data.frame(NULL)
  # iError_bounds <- c()
  for(id in nFP_estimate$ID) {
    # id = "N_E06"
    i_id <- which(levels(markedBounds_df$ID) == id)
    id_nFP <- nFP_estimate %>% filter(ID == id) %>% pull(nFP)
    id_singles <- singles_df %>% filter(ID == id)
    if(nrow(id_singles) > id_nFP) {
      i_error <- sample(1:nrow(id_singles), id_nFP, replace = F)
      id_singles$is_error[i_error] <- TRUE
    } else if(nrow(id_singles) > 0) {
      id_singles$is_error <- TRUE
    }
    singles_mrkErr <- rbind(singles_mrkErr, id_singles)
    # iSingles <- which(markedBounds_df$ID == id &
    #                     markedBounds_df$GT != "0/1" &
    #                     markedBounds_df$length == 1 &
    #                     markedBounds_df$LOH_cmplx == 0)
    # iLOH <- which(markedBounds_df$ID == id &
    #                 markedBounds_df$GT != "0/1")
    # if(remove_type == "singles"){
    #   nSingles <- length(iSingles)
    #   nError <- ifelse(nSingles >= id_nFP, id_nFP, nSingles)
    #   if(nSingles > 0) {
    #     iError_bounds_id <- sample(iSingles, nError)
    #   } else{
    #     iError_bounds_id <- c()
    #   }
    # 
    # } else {
    #   print("error")
    # }
    # else if(remove_type == "random") {
    #   iError_bounds_id <- c()
    #   nError <- 0
    #   nLOH <- sum(markedBounds_df$length[iLOH])
    #   nSample <- ifelse(nLOH >= id_nFP, id_nFP, nLOH)
    #   pLOH <-  markedBounds_df$length[iLOH]/nLOH
    #   for(s in 1:nSample) {
    #     iError_LOH <- sample(iLOH, 1, prob = pLOH)
    #     markedBounds_df$length[iError_LOH] <- markedBounds_df$length[iError_LOH] - 1
    #     iError_bounds_id <- c(iError_bounds_id, iError_LOH)
    #     if(markedBounds_df$length[iError_LOH] == 0) {
    #       iLOH <- iLOH[!iLOH == iError_LOH]
    #       pLOH <- pLOH[!iLOH == iError_LOH]
    #       nError <- nError + 1
    #     }
    #   }
    #   iError_bounds_id <- unique(iError_bounds_id)
    # }
    
    # iError_bounds <- c(iError_bounds, iError_bounds_id)
    # iError_bounds <- iError_bounds | iError_bounds_id
    
    # nErrorLOH_id <- data.frame(ID = id, nRmv = nError)
    # nErrorLOH <- rbind(nErrorLOH, nErrorLOH_id)
    if(i_id %% 10 == 0) {
      print(paste0(round(i_id/length(levels(markedBounds_df$ID)), 2)*100, "% done"))
    }
  }
  markedBounds_mrkErr_df <- rbind(singles_mrkErr, all_others_df)
  markedBounds_mrkErr_df <- markedBounds_mrkErr_df %>% arrange(ID, start_POSi)
  return(markedBounds_mrkErr_df)
}

MarkErrorLOH_pooled <- function(markedBounds_df, error_rate) {
  # markedBounds_df <- all_GT_bounds_merge
  # error_rate <- overall_F_Hom_rate
  rownames(markedBounds_df) <- NULL
  # Total number of markers sampled
  nSites_pool <- markedBounds_df %>% 
    summarise(nSites = sum(length))
  # Expected number of false homozygous calls
  nFP_estimate <- ceiling(nSites_pool$nSites * error_rate)
  
  markedBounds_df$is_error <- F
  # index LOHs, singletons, and doubletons
  i_LOH <- markedBounds_df$GT != "0/1"
  i_singles <- which(markedBounds_df$length == 1 & i_LOH)
  i_doubles <- which(markedBounds_df$length == 2 & i_LOH)
  
  
  if(length(i_singles) > nFP_estimate) {
    # if more singletons than false homs, randomly sample singletons as errors
    i_err_singles <- sample(i_singles, nFP_estimate, replace = F)
  } else if(length(i_singles) == nFP_estimate) {
    # if same number of singletons and false homs, all singletons are errors
    i_err_singles <- i_singles
  } else if(length(i_singles) < nFP_estimate &
            length(i_doubles) > nFP_doubles) {
    # if more false homs than singletons, all singletons are errors
    # and randomly sample doubletons for errors
      nFP_doubles <- nFP_estimate - length(i_singles)
    i_err_singles <- i_singles
    max_i_sampled <- 3
    while(max_i_sampled > 2) {
      # sample doublets with replacement. Only keep sampling if 
      # all doublets have been sampled <= 2 times
    i_err_doubles <- sample(i_doubles, nFP_doubles, replace = T)
    err_doubles_t <- table(i_err_doubles)
    max_i_sampled <- max(err_doubles_t)
    }
    is_single_err <- err_doubles_t == 1
    i_dbls_single_err <- as.numeric(names(err_doubles_t)[is_single_err])
    i_dbls_double_err <- as.numeric(names(err_doubles_t)[!is_single_err])
  } else {
    print("too many errors to accomodate")
  }
  markedBounds_df$length[i_dbls_single_err] <- 1
  markedBounds_df$is_error[c(i_err_singles, i_dbls_double_err)] <- T
  return(markedBounds_df)
}

MarkTerminalLOHs <- function(markedBounds_df, ancHet = LOH_SNPs) {
  # markedBounds_df <- all_LOHbounds_merge_EC
  if(class(markedBounds_df$ID) != "factor") {
    markedBounds_df$ID <- factor(markedBounds_df$ID)
  }
  markedBounds_df$isTerm <- F
  for (id in levels(markedBounds_df$ID)) {
    # id = "F_F06"
    i_id <- which(levels(markedBounds_df$ID) == id)
    id_bounds <- markedBounds_df %>% filter(ID == id)
    if (nrow(id_bounds) > 16) {
      id_chrBnds <- ancHet %>% 
        filter(ID == id) %>% group_by(CHROM) %>%
        summarise(chr_start = min(POSi, na.rm = T), chr_end = max(POSi, na.rm = T), .groups = "drop")
      id_chrMinMax <- id_bounds %>% 
        group_by(CHROM) %>% 
        summarize(minPOSi = min(start_POSi, na.rm = T), 
                  maxPOSi = max(end_POSi, na.rm = T), .groups = "drop")
      id_chrMinMax$minGT <- sapply(id_chrMinMax$minPOSi, 
                                   function(x) id_bounds$GT[id_bounds$start_POSi == x])
      id_chrMinMax$maxGT <- sapply(id_chrMinMax$maxPOSi, 
                                   function(x) id_bounds$GT[id_bounds$end_POSi == x])
      id_chrMin <- id_chrMinMax %>% filter(minGT != "0/1")
      id_chrMax <- id_chrMinMax %>% filter(maxGT != "0/1")
      startTerm <- sapply(id_chrMin$minPOSi, function(x) which(id_bounds$start_POSi == x))
      endTerm <- sapply(id_chrMax$maxPOSi, function(x) which(id_bounds$end_POSi == x))
      allTerm <- unlist(c(startTerm, endTerm))
      # allTerm <- c(unlist(ifelse(length(startTerm) > 0, startTerm, list(NULL))), 
      #              unlist(ifelse(length(endTerm) > 0, endTerm, list(NULL))))
      if(length(allTerm) > 0) {
        id_bounds$isTerm[allTerm] <- T
        cmplxTerm_k <- id_bounds %>% filter(isTerm == T & LOH_cmplx == 1) %>% pull(LOH_k)
        id_bounds$isTerm[id_bounds$LOH_k %in% cmplxTerm_k] <- T
      }
      markedBounds_df$isTerm[markedBounds_df$ID == id] <- id_bounds$isTerm
    }
    if(i_id %% 10 == 0) {
      print(paste0(round(i_id/length(levels(markedBounds_df$ID)), 2)*100, "% done"))
    }
  }
  return(markedBounds_df)
}

CountLOHevents <- function(markedBounds_df, omitError = F) {
  # markedBounds_df <- all_LOHbounds_merge_EC
  out_df <- data.frame(NULL)
  markedBounds_df$ID <- factor(markedBounds_df$ID)
  markedBounds_df$GT <- factor(markedBounds_df$GT, 
                               levels = c("0/0", "0/1", "1/1"))
  if (omitError == T) {
    markedBounds_df <- markedBounds_df %>% filter(is_error == F)
  }
  for (id in levels(markedBounds_df$ID)) {
    # id_df <- markedBounds_df %>% filter(ID == "F_A02")
    id_df <- subset(markedBounds_df, ID == id)
    n_LOHsmpl <- sum(id_df$GT != "0/1" & id_df$LOH_cmplx == 0)
    n_BY <- sum(id_df$GT == "0/0" & id_df$LOH_cmplx == 0)
    n_RM <- sum(id_df$GT == "1/1" & id_df$LOH_cmplx == 0)
    n_LOHcmplx <- max(id_df$LOH_k)
    id_LOHcount <- n_LOHsmpl + n_LOHcmplx
    n_BYcmplx <- 0
    n_RMcmplx <- 0
    n_mixCmplx <- 0
    if (max(id_df$LOH_k) > 0) {
      for (cx in 1:max(id_df$LOH_k)) {
        id_cmplx <- id_df %>% filter(LOH_k == cx)
        if (sum(id_cmplx$GT == "0/0") > 0 & sum(id_cmplx$GT == "1/1") == 0) {
          n_BYcmplx <- n_BYcmplx + 1
        }
        if (sum(id_cmplx$GT == "0/0") == 0 & sum(id_cmplx$GT == "1/1") > 0) {
          n_RMcmplx <- n_RMcmplx + 1
        }
        if (sum(id_cmplx$GT == "0/0") > 0 & sum(id_cmplx$GT == "1/1") > 0) {
          n_mixCmplx <- n_mixCmplx + 1
        }
      }
    }
    # nSmplTerm <- id_df %>% filter(isTerm == T, LOH_cmplx == 0) %>% nrow()
    # nCmplxTerm <- id_df %>% filter(isTerm == T, LOH_cmplx == 1) %>% distinct(LOH_k) %>% nrow()
    # nTerm <- nSmplTerm + nCmplxTerm
    # nInter <- id_LOHcount - nTerm
    # id_out <- data.frame(ID = id, n_LOH = id_LOHcount, 
    #                      n_smpl = n_LOHsmpl, n_BYsmpl = n_BY, n_RMsmpl = n_RM, 
    #                      n_cmplx = n_LOHcmplx, n_BYcmplx = n_BYcmplx, n_RMcmplx = n_RMcmplx, n_Mix = n_mixCmplx,
    #                      n_Term = nTerm, n_Inter = nInter)
    nTerm <- id_df %>% filter(isTerm == T) %>% nrow()
    nInter <- id_LOHcount - nTerm
    id_out <- data.frame(ID = id, n_LOH = id_LOHcount, 
                         n_smpl = n_LOHsmpl, n_BYsmpl = n_BY, n_RMsmpl = n_RM, 
                         n_cmplx = n_LOHcmplx, n_BYcmplx = n_BYcmplx, n_RMcmplx = n_RMcmplx, n_Mix = n_mixCmplx,
                         n_Term = nTerm, n_Inter = nInter)
    print(id)
    out_df <- rbind(out_df, id_out)
  }
  
  out_df$Tx <- factor(substr(out_df$ID, 1, 1))
  out_df$Tx_name <- Recode(out_df$Tx, "WT" = "N", "Cas9" = "H", "Drive" = "F")
  out_df$Rep <- factor(substr(out_df$ID, 4, 5))
  return(out_df)
}

find_indel_repeats <- function(seq) {
  # seq <- a_BY[2]
  if(is.na(seq) | nchar(seq) <= 3) {
    out <- F
  } else if(nchar(seq) <= 5) {
    out <- sum(rle(strsplit(seq, "(?<=.{1})", perl = TRUE)[[1]])$lengths >= 3) >= 1
  } else {
    out <- sum(c(rle(strsplit(seq, "(?<=.{1})", perl = TRUE)[[1]])$lengths >= 6, 
                 rle(strsplit(seq, "(?<=.{2})", perl = TRUE)[[1]])$lengths >= 3, 
                 rle(strsplit(substr(seq, 2, nchar(seq)), "(?<=.{2})", perl = TRUE)[[1]])$lengths >= 3,
                 rle(strsplit(seq, "(?<=.{3})", perl = TRUE)[[1]])$lengths >= 3, 
                 rle(strsplit(substr(seq, 2, nchar(seq)), "(?<=.{3})", perl = TRUE)[[1]])$lengths >= 3)) >= 1
  }
  return(out)
}

match_GT_alleles <- function(match_df = match_i_df, indels_df = indels_merge_same_n) {
  # Take indel df with split GT allele designations for each call set and 
  # columns for the final RM call alleles, as well as a map of the indicies 
  # of the alleles in the BY callset that match the RM callset
  for(posi in match_df$POSi){
    # posi = 1249076
    # posi = match_df$POSi[522]
    i_mis <- which(indels_df$POSi == posi)
    i_match <- which(match_df$POSi == posi)
    mis_1 <- as.numeric(indels_df$A_1_BYcall[i_mis])
    mis_1 <- as.numeric(indels_df$A_1_RMcall[i_mis])
    mis_2 <- as.numeric(indels_df$A_2_RMcall[i_mis])
    match_1 <- as.numeric(match_df[i_match, mis_1 + 2]) - 1
    match_2 <- as.numeric(match_df[i_match, mis_2 + 2]) - 1
    indels_df$A_1_final[i_mis] <- match_1
    indels_df$A_2_final[i_mis] <- match_2
  }
  return(indels_df)
}

perm_test <- function(df_in, cat_var, cat_names = NULL, response_var, test_stat = mean, 
                      n_perms = 10000, alpha = 0.05, alt_hyp = "two-tailed", rtrn = "p_df", include_matrix = F) {
  require(permute)
  # df_in = LOH_BPrate_SW_C
  # cat_var = "hs_C", cat_names = c("n", "y"), response_var = "p"
  # n_perms = 100, test_stat = mean
  
  # Ensure df_in is a data.frame
  df_in <- as.data.frame(df_in)

  # Ensure cat_var is a factor with no spurious levels
  df_in[, cat_var] <- factor(df_in[, cat_var])
  # df_in[, cat_var] <- droplevels(df_in[, cat_var])
  row.names(df_in) <- NULL
  # Create category name objects from "cat_names" or the first two levels of the "cat_var"
  if (length(cat_names) > 0 ) {
    cat1 <- cat_names[1]
    cat2 <- cat_names[2]
  } else {
    cat1 <- levels(df_in[, cat_var])[1]
    cat2 <- levels(df_in[, cat_var])[2]
  }
  # Indices of each category
  cat1_i <- df_in[, cat_var] == cat1
  cat2_i <- df_in[, cat_var] == cat2
  # remove levels not in cat_names
  df_in <- df_in[cat1_i | cat2_i, ]
  # Update indicies of each category
  cat1_i <- df_in[, cat_var] == cat1
  cat2_i <- df_in[, cat_var] == cat2
  # Create a vector of the categorical variables to index and shuffle
  cat_v <- as.vector(df_in[, cat_var])
  # Create empty premutation matrix
  dist_matrix <- matrix(nrow = n_perms, ncol = 2, dimnames = list(rows = NULL, cols = c("Trial", "Value")))
  
  if (is.function(test_stat)) {
    obs_val <- test_stat(df_in[cat2_i, response_var]) - test_stat(df_in[cat1_i, response_var])
  } 
  if (is.character(test_stat)) {
    if (test_stat == "chisq") {
      obs_val <- chisq.test(x = df_in[cat1_i, response_var], y = df_in[cat2_i, response_var])$statistic[[1]]
    } else if (test_stat == "chisq_hist") {
      count_1 <- df_in[cat1_i, response_var] %>% table(., dnn = list("n_rVar")) %>% as.data.frame(., responseName = "count")
      count_1$n_rVar <- levels(count_1$n_rVar)[count_1$n_rVar] %>% as.numeric()
      count_2 <- df_in[cat2_i, response_var] %>% table(., dnn = list("n_rVar")) %>% as.data.frame(., responseName = "count")
      count_2$n_rVar <- levels(count_2$n_rVar)[count_2$n_rVar] %>% as.numeric()
      n_rVar_seq <- list(min(c(count_1$n_rVar, count_2$n_rVar)):max(c(count_1$n_rVar, count_2$n_rVar)))
      names(n_rVar_seq) <- "n_rVar"
      count_1 <- merge(count_1, n_rVar_seq, by = "n_rVar", all = T, suffixes = c("1"))
      count_1_2 <- merge(count_1, count_2, by = "n_rVar", all = T, suffixes = c("1", "2"))
      count_1_2[is.na(count_1_2)] <- 0
      obs_val <- chisq.test(x = count_1_2$count1, y = count_1_2$count2)$statistic[[1]]
    } else if (test_stat == "t.test") {
      if (alt_hyp == "two-tailed"){
        alt = "two.sided"
      } else {
        alt = alt_hyp
      }
      obs_val <- t.test(x = df_in[cat1_i, response_var], y = df_in[cat2_i, response_var], 
                        alternative=alt, conf.level = 1-alpha)$statistic[[1]]
    } 
    else {
      print("test_stat not recognized")
      stop()
    }
  }

  # Perform permutation by shuffling categorical names among response values
  for (p in 1:n_perms) {
    # p = 1
    cat_v_shuff <- cat_v[shuffle(cat_v)]
    cat1_shuff_i <- cat_v_shuff == cat1
    cat2_shuff_i <- cat_v_shuff == cat2
    
    if (is.function(test_stat) == T) {
      shuff_val <- test_stat(df_in[cat1_shuff_i, response_var]) - test_stat(df_in[cat2_shuff_i, response_var])
    }
    if (is.character(test_stat)) {
      if (test_stat == "chisq") {
        shuff_val <- chisq.test(x = df_in[cat1_shuff_i, response_var], 
                                y = df_in[cat2_shuff_i, response_var])$statistic[[1]]
      }
      if (test_stat == "chisq_hist") {
        count_1 <- df_in[cat1_shuff_i, response_var] %>% table(., dnn = list("n_rVar")) %>% as.data.frame(., responseName = "count")
        count_1$n_rVar <- levels(count_1$n_rVar)[count_1$n_rVar] %>% as.numeric()
        count_2 <- df_in[cat2_shuff_i, response_var] %>% table(., dnn = list("n_rVar")) %>% as.data.frame(., responseName = "count")
        count_2$n_rVar <- levels(count_2$n_rVar)[count_2$n_rVar] %>% as.numeric()
        n_rVar_seq <- list(min(c(count_1$n_rVar, count_2$n_rVar)):max(c(count_1$n_rVar, count_2$n_rVar)))
        names(n_rVar_seq) <- "n_rVar"
        count_1 <- merge(count_1, n_rVar_seq, by = "n_rVar", all = T, suffixes = c("1"))
        count_1_2 <- merge(count_1, count_2, by = "n_rVar", all = T, suffixes = c("1", "2"))
        count_1_2[is.na(count_1_2)] <- 0
        shuff_val <- chisq.test(x = count_1_2$count1, y = count_1_2$count2)$statistic[[1]]
      }
      if (test_stat == "t.test") {
        shuff_val <- t.test(x = df_in[cat1_shuff_i, response_var], y = df_in[cat2_shuff_i, response_var], 
                            alternative=alt, conf.level = 1-alpha)$statistic[[1]]
      }
    }
    dist_matrix[p, 1] <- p
    dist_matrix[p, 2] <- shuff_val
  }
  # Calculate critical value from permuted data and alpha level according to chosen alternative
  # hypothesis
  if (alt_hyp == "two-tailed") {
    crit_val1 <- sort(dist_matrix[,2], decreasing = T)[floor(n_perms*alpha/2)]
    crit_val2 <- sort(dist_matrix[,2], decreasing = F)[floor(n_perms*alpha/2)]
    crit_val <- mean(abs(c(crit_val1, crit_val2)))
    reject_null <- abs(obs_val) > crit_val
    p_val <- 1 - sum(abs(obs_val) > abs(dist_matrix[,2]))/n_perms
  } 
  if (alt_hyp == "greater") {
    crit_val <- sort(dist_matrix[,2], decreasing = T)[floor(n_perms*alpha)]
    reject_null <- obs_val > crit_val
    p_val <- 1 - sum(obs_val > dist_matrix[,2])/n_perms
  }
  if (alt_hyp == "lesser") {
    crit_val <- sort(dist_matrix[,2], decreasing = F)[floor(n_perms*alpha)]
    reject_null <- obs_val < crit_val
    p_val <- 1 - sum(obs_val < dist_matrix[,2])/n_perms
  }
  # Calculate p-value as the proportion of the distribution as or more extreme than 
  # the observed data
  # p_val <- 1 - sum(abs(obs_val) > abs(dist_matrix[,2]))/n_perms
  
  # Place all values into list for outpu
  if(include_matrix) {
    results_list <- list(obsvStat = obs_val, 
                         critVal = crit_val, 
                         pVal = p_val,
                         rejectNull = reject_null, 
                         permMatrix = as.data.frame(dist_matrix))
  } else {
    results_list <- list(obsvStat = obs_val, 
                         critVal = crit_val, 
                         pVal = p_val,
                         rejectNull = reject_null)
  }
  
  if(rtrn == "p_df") {
      df_out <- data.frame(p_value = results_list[[3]])
      return(df_out)
    } else if(rtrn == "p") {
      return(results_list[[3]])
    } else {
    return(results_list)
  }
  
}

perm_test_DT <- function(DT, cat_var = NULL, response_var = NULL, test_stat = diff_means_DT, test_format = "factor",
                         n_perms = 1000, alpha = 0.05, alt_hyp = "two-tail") {
  # Requires a data.frame object with two columns, one with two factors for parsing the data
  # and one with the values to be compared
  if(length(cat_var) != 0 & length(response_var) != 0) {
    DT_cols <- c(cat_var, response_var)
    obs_DT <- DT[, ..DT_cols]
  } else {
    obs_DT <- DT[, 1:2]
  }
  names(obs_DT) <- c("cat_var", "rsp_var")
  obs_DT <- obs_DT[order(cat_var)]
  
  # Construct permutation data.table. Sample values from each group with replacement
  # and label with permutation number and individual number
  perms_AB <- obs_DT[, .(rsp_var = unlist(lapply(vector(mode = "list", length = n_perms),
                                                 function(x) .SD[sample(.N, .N, replace = F), rsp_var])),
                         cat_var = unlist(lapply(vector(mode = "list", length = n_perms), 
                                                 function(x) .SD[, cat_var])))][
                                                   , perm := rep(1:n_perms, each = .N/n_perms), by = "cat_var"][
                                                     , sample := rep(1:(.N), each = 1), by = c("cat_var", "perm")]
  
  # validate that permutations returning indepenedent samples of the pooled individuals for both groups
  # as.data.frame(perms_AB[perm %chin% 1:10]) %>% 
  #   ggplot() + geom_histogram(aes(x = rsp_var), binwidth = 1) + 
  #   facet_grid(cat_var~perm)
  
  # if function is from the stat tests, can feed in response and category variables into formula
  if(test_format == "formula") {
    # If test_stat is a {stats} function
    alt_hyp_map <- c("two.sided", "less", "greater")
    names(alt_hyp_map) <- c("two-tail", "left-tail", "right-tail")
    stat_alt <- alt_hyp_map[names(alt_hyp_map) == alt_hyp]
    obs_test_stat <- unlist(test_stat(rsp_var ~ cat_var, 
                                      data = obs_DT, 
                                      alternative = alt_hyp_map[alt_hyp][[1]])[
                                        c("statistic", "p.value")])
    
    test_stat_AB <- perms_AB[, lapply(.SD, function(x) test_stat(x ~ get(cat_var))$statistic), 
                             keyby = perm, .SDcols = response_var]
    setnames(test_stat_AB, response_var, "stat") # response variable becomes test statistic
  } else{
    obs_test_stat <- test_stat(obs_DT)
    test_stat_AB <- perms_AB[, .(stat = test_stat(.SD)), 
                             keyby = perm, .SDcols = c("cat_var", "rsp_var")][order(stat)]
  }
  if(alt_hyp == "two-tail") {
    i_test_low <- floor(n_perms * (alpha/2))
    i_test_hi <- ceiling(n_perms * (1 - alpha/2))
    crit_low <- test_stat_AB[i_test_low, stat]
    crit_hi <- test_stat_AB[i_test_hi, stat]
    crit_val <- mean(abs(c(crit_low, crit_hi)))
    sig <- abs(obs_test_stat[[1]]) > crit_val
    p_val <- 1 - sum(abs(obs_test_stat[[1]]) > test_stat_AB[order(abs(stat))][, abs(stat)])/n_perms
    
  } else if(alt_hyp == "left-tail") {
    i_test_low <- floor(n_perms * alpha)
    crit_low <- test_stat_AB[i_test_low, stat]
    sig <- obs_test_stat[[1]] < crit_low
    p_val <- test_stat_AB[stat < obs_test_stat[[1]]][, .N]/n_perms
  } else if(alt_hyp == "right-tail") {
    i_test_hi <- ceiling(n_perms * (1 - alpha))
    crit_hi <- test_stat_AB[i_test_hi, stat]
    sig <- obs_test_stat[[1]] > crit_hi
    p_val <- test_stat_AB[stat > obs_test_stat[[1]]][, .N]/n_perms
  } else {
    print("alt_hyp must be one of c(two-tail, left-tail, right-tail))")
  }
  return(p_val)
}

diff_means_DT <- function(DT_in) {
  # takes a data.table with categorical and response variable columns
  # and takes the difference in means 
  DT_in[, .(m = mean(rsp_var)), by = cat_var][, diff(m)]
} 


BHcorrection <- function(p_df, p_col = "pVal", FDR = 0.05) {
  # p_df <- rstest_df
  rownames(p_df) <- NULL
  p_df$i_original <- 1:nrow(p_df)
  p_df_order <- p_df %>% arrange(across(contains(p_col)))
  n_tests <- nrow(p_df_order)
  p_df_order$i <- 1:n_tests
  p_df_order$p_crit_BH <- FDR * (p_df_order$i/n_tests)
  iCritical <- p_df_order[, p_col] < p_df_order$p_crit_BH
  p_df_order$rejectNull_BH <- 0
  if(sum(iCritical) > 0) {
    critical_row <- max(as.numeric(rownames(p_df_order[iCritical, ])))
    p_df_order$rejectNull_BH[1:critical_row] <- 1
  }
  p_df_out <- p_df_order %>% arrange(i_original) %>% select(!c(i_original, i))
  return(p_df_out)
}

get_LOH_coordinates <- function(POSi_df, Het = F, tLOH_offset = 10000) {
  # POSi_df <- all_LOHbounds_merge_NS
  POSi_df$est_start_POS <- ConvertPosIndicies(POSi_df, pos_col = "est_start", index_out = "POS")
  POSi_df$est_end_POS <- ConvertPosIndicies(POSi_df, pos_col = "est_end", index_out = "POS")
  
  if(!Het) {
    POSi_df$dist_cent_start <- ChromosomeCoordinates(POSi_df, POSi_col = "est_start") %>% 
      pull(dist_cent)
    POSi_df$dist_cent_end <- ChromosomeCoordinates(POSi_df, POSi_col = "est_end") %>% 
      pull(dist_cent)
    
    POSi_df <- POSi_df %>% mutate(est_mid_POS = (est_start_POS + est_end_POS)/2)
    
    POSi_df$ID_CHROM <- paste0(POSi_df$ID, "_", POSi_df$CHROM)
    
    i_Term <- POSi_df$isTerm
    i_large_term <- i_Term & POSi_df$est_length >= 2 * tLOH_offset
    head_term <- POSi_df$dist_cent_start < 0 & POSi_df$dist_cent_end < 0
    # tail_term <- POSi_df$dist_cent_start > 0 & POSi_df$dist_cent_end > 0
    cross_h_term <- POSi_df$dist_cent_start < 0 & POSi_df$dist_cent_end > 0
    cross_ID_CHROM <- POSi_df[i_Term & cross_h_term, ] %>% pull(ID_CHROM)
    # cross_ID_CHROM <- POSi_df[i_Term & cross_h_term & !i_wh_chrom, ] %>% pull(ID_CHROM)
    
    cross_h_term_ID <- POSi_df %>%
      ungroup() %>%
      filter(i_large_term) %>%
      # mutate(ID_CHROM = paste0(ID, "_", CHROM)) %>% 
      filter(ID_CHROM %in% cross_ID_CHROM) %>% distinct(ID_CHROM, .keep_all = T) %>% 
      summarize(ID_CHROM = ID_CHROM, head_term = ifelse(est_start_POS == 1, T, F))
    
    # Correct POS for terminal LOH that do not cross the centromere
    POSi_df$est_mid_POS[i_large_term] <- ifelse(head_term[i_large_term], 
                                                POSi_df$est_end_POS[i_large_term] - tLOH_offset, 
                                                POSi_df$est_start_POS[i_large_term] + tLOH_offset)
    
    # Correct POS for terminal LOH that do cross the centromere
    i_crct <- i_large_term & cross_h_term
    # i_crct <- i_Term & cross_h_term & !i_wh_chrom
    
    POSi_df$est_mid_POS[i_crct][cross_h_term_ID$head_term] <- 
      POSi_df$est_end_POS[i_crct][cross_h_term_ID$head_term] - tLOH_offset
    POSi_df$est_mid_POS[i_crct][!cross_h_term_ID$head_term] <- 
      POSi_df$est_start_POS[i_crct][!cross_h_term_ID$head_term] + tLOH_offset
    POSi_df$est_mid <- ConvertPosIndicies(POSi_df, pos_col = "est_mid_POS", index_out = "POSi")
    POSi_df$dist_cent_mid <- ChromosomeCoordinates(POSi_df, POSi_col = "est_mid") %>% 
      pull(dist_cent)
    
    POSi_df$dist_term <- 0
    i_cent_neg <- POSi_df$dist_cent_mid < 0
    POSi_df$dist_term[i_cent_neg] <- POSi_df$est_mid_POS[i_cent_neg]
    POSi_df$dist_term[!i_cent_neg] <- operate_by_factor_match(chrom_bound_BY[, c("CHROM", "End")],
                                                              POSi_df[!i_cent_neg, c("CHROM", "est_mid")],
                                                              .fun = function(x, y) x - y)
    POSi_df$fract_dist <- 0
    POSi_df$fract_dist[i_cent_neg] <- operate_by_factor_match(chrom_arms[, c("CHROM", "arm_1_cent")],
                                                              POSi_df[i_cent_neg, c("CHROM", "dist_cent_mid")],
                                                              .fun = function(x, y) y/x)  
    POSi_df$fract_dist[!i_cent_neg] <- operate_by_factor_match(chrom_arms[, c("CHROM", "arm_2_cent")],
                                                               POSi_df[!i_cent_neg, c("CHROM", "dist_cent_mid")],
                                                               .fun = function(x, y) y/x)  
    
    POSi_df$dist_cent_mid_abs <- POSi_df$dist_cent_mid %>% round() %>% abs()
    POSi_df$type <- factor(ifelse(POSi_df$isTerm, "Terminal", "Interstitial"))
  }
  return(POSi_df)
}


SliderCalc <- function(df, data_col, index_col, factor_col, 
                       window_size, slide_interval, chrom_win = F, summary_stat, ...) {
  ## Takes in a dataframe containing a column of data values to be summarized within a window, an optional
  ## column of alternate indicies to apply the sliding window to, and a statistic to calculate within each
  ## window. Outputs a new dataframe with the summarized values and window indicies. The "summary_stat" can
  ## be one of the base R stats such as mean(), median(), min(), ect. which operates on a vector and outputs
  ## a single value
  # df <- POSi_data_wHet
  df <- as.data.frame(df)

  if(missing("slide_interval")) {
    slide_interval <- window_size
  }
  if(missing("index_col")) {
    pos_min <- 1
    pos_max <- nrow(df)
  } else if(index_col == "POS") {
    chrom_win <- T
  } else {
    pos_min <- min(df[, index_col])
    pos_max <- max(df[, index_col])
  }
   
  # Create windows. If chrom_win = T, makes windows according to chromosome 
  # boundaries and window_size
  # Otherwise, makes windows according to min and max positions and 
  # window_size
  if(class(window_size) == "data.frame") {
    i_s <- grep("start", colnames(window_size), ignore.case = T)
    i_e <- grep("end", colnames(window_size), ignore.case = T)
    i_c <- grep("chr", colnames(window_size), ignore.case = T)
    start_index <- window_size[, i_s]
    end_index <- window_size[, i_e]
    if(length(i_c) != 0){
      chr_name <-  window_size[, i_c]
    }
  } else {
    if(sum(pos_max < window_size | pos_max < slide_interval) != 0) {
      return("Window too big for data")
    }
    if(chrom_win) {
      df$CHROM <- factor(df$CHROM)
      # gen_wndws_df <- data.frame(NULL)
      start_index <- c()
      end_index <- c()
      chr_index <- c()
      if(index_col == "POS") {
        for(ch in 1:nrow(chrom_bound_BY)) {
          # ch = 1
          chrom_wndws_start <- c(seq(from = 1, 
                                     to = chrom_bound_BY$length[ch], 
                                     by = slide_interval))
          chrom_wndws_end <- chrom_wndws_start + window_size - 1
          make_up <- (max(chrom_wndws_end) - chrom_df$length[ch])/(length(chrom_wndws_end) - 1)
          offsets <- round(cumsum(rep(make_up, (length(chrom_wndws_end) - 1))))
          chrom_wndws_start[-1] <- chrom_wndws_start[-1] - offsets
          chrom_wndws_end[-1] <- chrom_wndws_end[-1] - offsets
          chrom_chr <- rep(ch, length(chrom_wndws_start))
          start_index <- c(start_index, chrom_wndws_start)
          end_index <- c(end_index, chrom_wndws_end)
          chr_index <- c(chr_index, chrom_chr)
        }
        chr_name <- str_pad(chr_index, width = 2, pad = "0")
        chr_window_df <- data.frame(CHROM = chr_name, 
                                    start = start_index, end = end_index, 
                                    chr_index = chr_index)
      } else if(index_col == "POSi") {
        for(ch in 1:nrow(chrom_bound_BY)) {
          # ch = 16
          chrom_wndws_start <- c(seq(from = chrom_df$Start[ch], 
                                     to = chrom_df$End[ch], 
                                     by = slide_interval))
          chrom_wndws_end <- chrom_wndws_start + window_size - 1
          make_up <- (max(chrom_wndws_end) - chrom_df$End[ch])/(length(chrom_wndws_end) - 1)
          offsets <- round(cumsum(rep(make_up, (length(chrom_wndws_end) - 1))))
          chrom_wndws_start[-1] <- chrom_wndws_start[-1] - offsets
          chrom_wndws_end[-1] <- chrom_wndws_end[-1] - offsets
          chrom_chr <- rep(ch, length(chrom_wndws_start))
          start_index <- c(start_index, chrom_wndws_start)
          end_index <- c(end_index, chrom_wndws_end)
          chr_index <- c(chr_index, chrom_chr)
      }
      chr_name <- str_pad(chr_index, width = 2, pad = "0")
      chr_window_df <- data.frame(CHROM = chr_name, 
                                  start = start_index, end = end_index, 
                                  chr_index = chr_index)
      } else {
        print("Chromosome windows must index POS or POSi")
      }
    } else {
      start_index <- seq(from = pos_min, to = pos_max, by = slide_interval)
      end_index <- start_index + window_size - 1
    } 
  }
  if (missing("factor_col")) {
    if(missing("index_col")) {
      pos_idx <- 1:nrow(df)
      win_smry <- c()
      n_vec <- c()
      for (w in seq_along(start_index)) {
        # w = 18
        win_i <- pos_idx >= start_index[w] & pos_idx <= end_index[w]
        n_ele <- sum(win_i)
        win_smry <- c(win_smry, summary_stat(df[win_i, data_col], na.rm = T))
        n_vec <- c(n_vec, n_ele)
      }
      slider_out <- data.frame(start = start_index, end = end_index, 
                               n_elements = n_vec, value = win_smry)
      colnames(slider_out)[4] <- data_col
    } else if(index_col == "POS") {
      win_smry <- c()
      n_vec <- c()
      for (w in 1:nrow(chr_window_df)) {
        # w = 1
        win_i <- df[, index_col] >= chr_window_df$start[w] & 
          df[, index_col] <= chr_window_df$end[w] & 
          df$CHROM == chr_window_df$CHROM[w]
        n_ele <- sum(win_i)
        win_smry <- c(win_smry, summary_stat(df[win_i, data_col], na.rm = T))
        n_vec <- c(n_vec, n_ele)
      }
      slider_out <- cbind(chr_window_df[, -4], n_elements = n_vec, value = win_smry)
      colnames(slider_out)[5] <- data_col
    } else {
      pos_idx <- df[, index_col]
      win_smry <- c()
      n_vec <- c()
      for (w in seq_along(start_index)) {
        # w = 18
        win_i <- pos_idx >= start_index[w] & pos_idx <= end_index[w]
        n_ele <- sum(win_i)
        win_smry <- c(win_smry, summary_stat(df[win_i, data_col], na.rm = T))
        n_vec <- c(n_vec, n_ele)
      }
      if(chrom_win) {
        slider_out <- data.frame(start = start_index, end = end_index, 
                                 CHROM = chr_name,
                                 n_elements = n_vec, value = win_smry)
        colnames(slider_out)[5] <- data_col
      } else {
        slider_out <- data.frame(start = start_index, end = end_index, 
                                 n_elements = n_vec, value = win_smry)
        colnames(slider_out)[4] <- data_col
      }
    }
    return(slider_out)
  }
  else {
    df[, factor_col] <- droplevels(df[, factor_col])
    slider_out <- data.frame(NULL)
    for (f in levels(df[,factor_col])) {
      df_sub <- subset(df, df[,factor_col] == f)
      if(missing("index_col")) {
        pos_idx <- 1:nrow(df_sub)
        win_smry <- c()
        n_vec <- c()
        for (w in seq_along(start_index)) {
          # w = 18
          win_i <- pos_idx >= start_index[w] & pos_idx <= end_index[w]
          n_ele <- sum(win_i)
          win_smry <- c(win_smry, summary_stat(df_sub[win_i, data_col], na.rm = T))
          n_vec <- c(n_vec, n_ele)
        }
        smry_df <- data.frame(factor_col = f, start = start_index, end = end_index, 
                              n_elements = n_vec, value = win_smry)
        colnames(smry_df)[c(1, 5)] <- c(factor_col, data_col)
      } else if(index_col == "POS") {
        win_smry <- c()
        n_vec <- c()
        for (w in 1:nrow(chr_window_df)) {
          # w = 1
          win_i <- df[, index_col] >= chr_window_df$start[w] & 
            df[, index_col] <= chr_window_df$end[w] & 
            df$CHROM == chr_window_df$CHROM[w]
            n_ele <- sum(win_i)
            win_smry <- c(win_smry, summary_stat(df[win_i, data_col], na.rm = T))
            n_vec <- c(n_vec, n_ele)
        }
        smry_df <- cbind(factor_col = f, chr_window_df[, -4], n_elements = n_vec, value = win_smry)
        colnames(smry_df)[c(1, 6)] <- c(factor_col, data_col)
      } else {
        pos_idx <- df_sub[, index_col]
        win_smry <- c()
        n_vec <- c()
        for (w in seq_along(start_index)) {
          # w = 18
          win_i <- pos_idx >= start_index[w] & pos_idx <= end_index[w]
          n_ele <- sum(win_i)
          win_smry <- c(win_smry, summary_stat(df_sub[win_i, data_col], na.rm = T))
          n_vec <- c(n_vec, n_ele)
        }
        if(chrom_win) {
          smry_df <- data.frame(factor_col = f, start = start_index, end = end_index, 
                                CHROM = chr_name,
                                n_elements = n_vec, value = win_smry)
          colnames(smry_df)[c(1, 6)] <- c(factor_col, data_col)
          
        } else {
          smry_df <- data.frame(factor_col = f, start = start_index, end = end_index, 
                                n_elements = n_vec, value = win_smry)
          colnames(smry_df)[c(1, 5)] <- c(factor_col, data_col)
        }
      }
      slider_out <- rbind(slider_out, smry_df)
    }
    if(chrom_win) {
      slider_out$CHROM <- factor(slider_out$CHROM)
    }
    return(slider_out)
  }
}

GetRateAndData <- function(ancHet_in, by_GT = F, include_het = T, three_sites = F) {
  require(pammtools)
  # Takes the "all_ancHet" dataframe and calculates the LOH state and data presence across
  # every position in the genome. t_i representst the state of each site
  # r_i can have values 0 (het) or 1, 3 (hom)
  # d_i can have values 0 (no data) or 1 (data)
  # t_i is a combination of r_i and d_i
  # id = "N_G01"
  all_POSi <- ancHet_in %>% select(POSi, CHROM, POS) %>% 
    distinct() %>% arrange(POSi)
  t_runs_all <- data.frame(NULL)
  for (id in unique(ancHet_in$ID)) {
    # id = "F_G12"
    # ancHet_in = all_ancHet
    id_RMxBY <- ancHet_in %>% filter(ID == id) %>% 
      select(ID, POSi, CHROM, POS, GT)
    # Merge clone sites with master list of all sites with data across clones
    id_RMxBY <- merge(all_POSi, id_RMxBY, 
                      by = c("POSi", "CHROM", "POS"), 
                      sort = T, all = T) 
    if(by_GT == F) {
      id_RMxBY <- id_RMxBY %>% mutate(r_i = as.numeric(GT != "het"), d_i = as.numeric(!is.na(GT)))
      # Data (t_i) can be in three states
      # 0 indicates no data
      # 1 indicates LOH site
      # 2 indicates het site
      id_RMxBY$t_i <- 2
      id_RMxBY$t_i[id_RMxBY$r_i == 1 & id_RMxBY$d_i == 1] <- 1
      id_RMxBY$t_i[id_RMxBY$d_i == 0] <- 0
    }
    if (by_GT == T) {
      id_RMxBY <- id_RMxBY %>% mutate(r_i = as.numeric(GT), d_i = as.numeric(!is.na(GT)))
      # Data can be in four states
      # 0 indicates no data
      # 1 indicates BY site
      # 2 indicates het site
      # 3 indicates RM site
      id_RMxBY$t_i <- 2
      id_RMxBY$t_i[id_RMxBY$r_i == 1 & id_RMxBY$d_i == 1] <- 1
      id_RMxBY$t_i[id_RMxBY$r_i == 3 & id_RMxBY$d_i == 1] <- 3
      id_RMxBY$t_i[id_RMxBY$d_i == 0] <- 0
    }
    # Detect runs of t_i and create a dataframe of positions
    t_runs_df <- data.frame(NULL)
    for (ch in levels(id_RMxBY$CHROM)) {
      # ch = "12"
      t_runs <- rle(id_RMxBY$t_i[id_RMxBY$CHROM == ch])
      t_runs_ch <- data.frame(value = t_runs$values, 
                              length = t_runs$lengths,
                              CHROM = ch)
      t_runs_df <- rbind(t_runs_df, t_runs_ch)
    }
    t_runs_df$start <- cumsum(t_runs_df$length) - t_runs_df$length + 1
    t_runs_df$end <-  cumsum(t_runs_df$length)
    t_runs_df$start_POSi <- id_RMxBY$POSi[t_runs_df$start]
    t_runs_df$end_POSi <- id_RMxBY$POSi[t_runs_df$end]
    # Estimate boundaries of type events
    t_runs_df$est_start <- 0
    t_runs_df$est_end <- 0
    if (three_sites == T) {
      t_runs_df$short <- 0
      t_runs_df$short[t_runs_df$value %in% c(1,3) & t_runs_df$length < 3] <- 1
      t_runs_df <- t_runs_df %>% filter(short == 0) %>% select(-short)
    }
    for (i in 1:nrow(t_runs_df)) {
      # i = 2
      obsv_start <- t_runs_df$start_POSi[i]
      obsv_end <- t_runs_df$end_POSi[i]
      if (i == 1) {
        t_runs_df$est_start[i] <- obsv_start
      } else {
        t_runs_df$est_start[i] <- ceiling((obsv_start + t_runs_df$end_POSi[i - 1]) / 2)
      } 
      if (i == nrow(t_runs_df)) {
        t_runs_df$est_end[i] <- obsv_end
      } else {
        t_runs_df$est_end[i] <- floor((obsv_end + t_runs_df$start_POSi[i + 1] - 1) / 2)
      }
    }
    t_runs_df$est_length <- t_runs_df$est_end - t_runs_df$est_start + 1
    if(include_het == T) {
      t_runs_df$ID <- id
      t_runs_all <- rbind(t_runs_all, t_runs_df)
    } else {
      t_pos_df <- t_runs_df %>% filter(value != 2)
      t_pos_df$ID <- id
      t_runs_all <- rbind(t_runs_all, t_pos_df)
    }
    print(id)
  }
  t_runs_all <- CategoriesFromID(t_runs_all)
  return(t_runs_all)
}

LOHrateSW <- function(t_pos_df, by_GT = F, win = 50000, slide_dist = 50000, 
                      chrom_df = chrom_bound_BY, rm_zero_valids = F) {
  # t_pos_df <- POSi_data_wHet
  if(sum(t_pos_df$value == 2) == 0) {
    print("Heterozygous markers required for analysis")
    stop()
  }
  if (is.numeric(win)) {
    gen_wndws_df <- data.frame(NULL)
    for(ch in 1:nrow(chrom_df)) {
      # ch = 1
      chrom_wndws_start <- c(seq(from = chrom_df$Start[ch], 
                           to = chrom_df$End[ch], 
                           by = slide_dist))
      chrom_wndws_end <- chrom_wndws_start + win - 1
      make_up <- (max(chrom_wndws_end) - chrom_df$End[ch])/(length(chrom_wndws_end) - 1)
      offsets <- round(cumsum(rep(make_up, (length(chrom_wndws_end) - 1))))
      chrom_wndws_start[-1] <- chrom_wndws_start[-1] - offsets
      chrom_wndws_end[-1] <- chrom_wndws_end[-1] - offsets
      chrom_wndws <- data.frame(start_POSi = chrom_wndws_start, end_POSi = chrom_wndws_end)
      # if (ch == nrow(chrom_df)) {
      #   chrom_wndws <- c(chrom_wndws, chrom_df[ch, 2])
      # }
      chrom_wndws_df <- data.frame(start_POSi = chrom_wndws$start_POSi,
                                   end_POSi = chrom_wndws$end_POSi,
                                   CHROM = as.character(chrom_df$CHROM[ch]),
                                   start_POS = chrom_wndws$start_POSi - chrom_df$Start[ch] + 1,
                                   end_POS = chrom_wndws$end_POSi - chrom_df$Start[ch] + 1)
      gen_wndws_df <- rbind(gen_wndws_df, chrom_wndws_df)
    }
  } else {
    if (win == "CHROM") {
      # win = "CHROM"
      gen_wndws_df <- data.frame(start_POSi = chrom_df[, "Start"],
                                 end_POSi = chrom_df[, "End"],
                                 CHROM = as.character(chrom_df$CHROM),
                                 start_POS = 1,
                                 end_POS = chrom_df[, "length"],)
      
    } else if (win == "CHROM/3") {
      gen_wndws_df <- data.frame(NULL)
      ch_len <- chrom_df$length
      for(ch in 1:16) {
        # ch = 2
        chrom_wndws_start <- c(chrom_df[ch, 1],
                         chrom_df[ch, 1] + ceiling(ch_len[ch]/3), 
                         chrom_df[ch, 1] + ceiling(ch_len[ch]*2/3))
        chrom_wndws_end <- c((chrom_wndws_start[-1] - 1),chrom_df$End[ch])
        chrom_wndws <- data.frame(start_POSi = chrom_wndws_start, end_POSi = chrom_wndws_end)
        chrom_wndws_df <- data.frame(start_POSi = chrom_wndws$start_POSi,
                                     end_POSi = chrom_wndws$end_POSi,
                                     CHROM = as.character(chrom_df$CHROM[ch]), 
                                     CHROM3 = paste0(as.character(chrom_df$CHROM[ch]), "_", 1:3),
                                     start_POS = chrom_wndws$start_POSi - chrom_df$Start[ch] + 1,
                                     end_POS = chrom_wndws$end_POSi - chrom_df$Start[ch] + 1)
        gen_wndws_df <- rbind(gen_wndws_df, chrom_wndws_df)
      }
      gen_wndws_df <- rbind(gen_wndws_df, chrom_wndws_df)
    } else {
      print("Window not recognized")
    }
  }
  # t_pos_df <- all_LOHbounds_noGT %>% filter(! (value == 1 & isTerm == F))
  full_sw_df <- data.frame(NULL)
  t_pos_df$ID <- factor(t_pos_df$ID)
  for (id in levels(t_pos_df$ID)) {
    # id = "F_A05"
    # t_pos_id <- POSi_data_wHet %>% filter(ID == "F_A07")
    t_pos_id <- t_pos_df %>% filter(ID == id)
    r_het_df <- t_pos_id %>% filter(value == 2)
    r_pos_df <- t_pos_id %>% filter(value == 1 | value == 3)
    # r_neg_df <- t_pos_id %>% filter(value == 0)
    # id_df_o <- id_df
    if(win == "CHROM/3") {
      id_df <- data.frame(ID = id, 
                          start_POSi = 0, end_POSi = 0,
                          CHROM = "C", CHROM3 = "C3", 
                          start_POS = 0, end_POS = 0,
                          n_sites = 0, LOH_rate = 0, prop_valid = 0)
    } else {
      id_df <- data.frame(ID = id, 
                          start_POSi = 0, end_POSi = 0,
                          CHROM = "C", 
                          start_POS = 0, end_POS = 0,
                          n_sites = 0, LOH_rate = 0, prop_valid = 0)
    }
    for(w in 1:(nrow(gen_wndws_df))) {
      # w = nrow(gen_wndws_df)-1
      # w <- 5
      # Calculate the number of t_i = 1 (Hom) sites for each window
      r_out_idx <- r_pos_df$est_start <= gen_wndws_df$start_POSi[w] & 
        r_pos_df$est_end >= gen_wndws_df$end_POSi[w]
      if (sum(r_out_idx) == 0) {
        r_in_idx <- r_pos_df$est_start >= gen_wndws_df$start_POSi[w] & 
          r_pos_df$est_end <= gen_wndws_df$end_POSi[w]
        n_r_in <- sum(r_pos_df$est_length[r_in_idx])
        r_lead_idx <- r_pos_df$est_start < gen_wndws_df$start_POSi[w] & 
          r_pos_df$est_end >= gen_wndws_df$start_POSi[w] &
          r_pos_df$est_end <= gen_wndws_df$end_POSi[w]
        if(sum(r_lead_idx) > 0) {
          n_r_lead <- r_pos_df$est_end[r_lead_idx] - gen_wndws_df$start_POSi[w] + 1
        } else {
          n_r_lead <- 0
        }
        r_trail_idx <- r_pos_df$est_start >= gen_wndws_df$start_POSi[w] & 
          r_pos_df$est_start < gen_wndws_df$end_POSi[w] & 
          r_pos_df$est_end > gen_wndws_df$end_POSi[w]
        if(sum(r_trail_idx) > 0) {
          n_r_trail <- gen_wndws_df$end_POSi[w] - r_pos_df$est_start[r_trail_idx]
        } else {
          n_r_trail <- 0
        }
        n_r_wndw <- sum(n_r_in, n_r_lead, n_r_trail)
      } else {
        n_r_wndw <- gen_wndws_df$end_POSi[w] - gen_wndws_df$start_POSi[w]
      }
      # w = 249
      # Calculate the number of t_i = 2 (Het) sites for each window
      neg_out_idx <- r_het_df$est_start <= gen_wndws_df$start_POSi[w] & r_het_df$est_end >= gen_wndws_df$end_POSi[w]
      if (sum(neg_out_idx) == 0) {
        neg_in_idx <- r_het_df$est_start >= gen_wndws_df$start_POSi[w] & r_het_df$est_end <= gen_wndws_df$end_POSi[w]
        n_neg_in <- sum(r_het_df$est_length[neg_in_idx])
        neg_lead_idx <- r_het_df$est_start <= gen_wndws_df$start_POSi[w] & 
          r_het_df$est_end >= gen_wndws_df$start_POSi[w] &
          r_het_df$est_end <= gen_wndws_df$end_POSi[w]
        n_neg_lead <- r_het_df$est_end[neg_lead_idx] - gen_wndws_df$start_POSi[w] + 1
        neg_trail_idx <- r_het_df$est_start >= gen_wndws_df$start_POSi[w] & 
          r_het_df$est_start <= gen_wndws_df$end_POSi[w] & 
          r_het_df$est_end >= gen_wndws_df$end_POSi[w]
        n_neg_trail <- gen_wndws_df$end_POSi[w] - r_het_df$est_start[neg_trail_idx]
        n_neg_wndw <- sum(n_neg_in, n_neg_lead, n_neg_trail)
      } else {
        n_neg_wndw <- gen_wndws_df$end_POSi[w] - gen_wndws_df$start_POSi[w]
      }
      n_site <- n_r_wndw + n_neg_wndw
      LOHrate <- n_r_wndw / n_site
      prop_d <- n_site / (gen_wndws_df$end_POSi[w] - gen_wndws_df$start_POSi[w])
      if(win == "CHROM/3") {
        id_df[w, ] <- list(id, gen_wndws_df$start_POSi[w], gen_wndws_df$end_POSi[w],
                           gen_wndws_df$CHROM[w], gen_wndws_df$CHROM3[w], 
                           gen_wndws_df$start_POS[w], gen_wndws_df$end_POS[w], 
                           n_site, LOHrate, prop_d)
      } else {
        id_df[w, ] <- list(id, gen_wndws_df$start_POSi[w], gen_wndws_df$end_POSi[w],
                           gen_wndws_df$CHROM[w], 
                           gen_wndws_df$start_POS[w], gen_wndws_df$end_POS[w], 
                           n_site, LOHrate, prop_d)
      }
    }
    if(rm_zero_valids) {
      id_df <- id_df %>% filter(prop_valid > 0)
    }
    
    full_sw_df <- rbind(full_sw_df, id_df)
    print(id)
  }
  full_sw_df <- CategoriesFromID(full_sw_df)
  return(full_sw_df)
}

se <- function(x, ...) {
  sqrt(var(x[!is.na(x)], ...)/length(x[!is.na(x)]))
}

rateMeans_xTx <- function(rates_df, per_gen = F, gens = n_gens, min_n = F, min_prop_valid = 0.25) {
  # rates_df <- all_LOHrates_SW50k
  all_LOHrates_sw <- data.frame(NULL)
  rates_df <- rates_df %>% ungroup() %>% dplyr::arrange(Tx_name, POSi)
  if(is.numeric(min_prop_valid)) {
    rates_df <- rates_df %>% filter(prop_valid > min_prop_valid)
  }
  for (tx in levels(rates_df$Tx_name)) {
    # tx <- "WT"
    LOHrates_wide <- rates_df %>% 
      filter(Tx_name == tx) %>% 
      select(ID, POSi, CHROM, POS, LOH_rate, n_sites) %>% 
      pivot_wider(names_from = ID, values_from = c(LOH_rate, n_sites))
    
    LOHrates_n_valid <- apply(LOHrates_wide %>% 
                                select(starts_with("LOH_rate")), 
                              MARGIN = 1, function(x) sum(!is.na(x)))
    if(is.numeric(min_n)) {
    LOHrates_wide <- LOHrates_wide[LOHrates_n_valid >= min_n, ]
    }
    LOHmean_v <- LOHrates_wide %>% 
      select(starts_with("LOH_rate")) %>% 
      apply(., mean, MARGIN = 1, na.rm = T)
    LOHse_v <- LOHrates_wide %>% 
      select(starts_with("LOH_rate")) %>% 
      apply(., se, MARGIN = 1, na.rm = T)
    if(per_gen) {
      LOHmean_v <- LOHmean_v/gens
      LOHse_v <- LOHse_v/gens
    }
    meanSites_v <- LOHrates_wide %>% 
      select(starts_with("n_sites")) %>% 
      apply(., mean, MARGIN = 1, na.rm = T)
    LOHrate_SW <- cbind(Tx_name = tx, 
                        LOHrates_wide[,c(1:3)], 
                        rate = LOHmean_v,
                        std_err = LOHse_v,
                        n_sites = meanSites_v,
                        n_clones = LOHrates_n_valid)
    
    all_LOHrates_sw <- rbind(all_LOHrates_sw, LOHrate_SW)
  }
  all_LOHrates_sw$Tx_name <- factor(all_LOHrates_sw$Tx_name, levels = c("WT", "Cas9", "Drive"))
  all_LOHrates_sw <- all_LOHrates_sw %>% group_by(Tx_name) %>% mutate(end_POSi = lead(POSi) - 1)
  all_LOHrates_sw$end_POSi[is.na(all_LOHrates_sw$end_POSi)] <- g_length
  all_LOHrates_sw$CHROM <- factor(all_LOHrates_sw$CHROM)
  return(all_LOHrates_sw)
}

fill_zeros <- function(count_df, group_col = "Tx_name", factor_col = "n_SNMs", fill_col = "n_muts") {
  count_df <- count_df %>% as.data.frame()
  max_n <- max(count_df[, factor_col])
  ns <- 0:max_n
  # max_ns <- count_df %>% group_by_(.dots = list(group_col)) %>% summarise(across(all_of(factor_col), max))
  for(tx in unique(count_df[, group_col])) {
    # tx <- "Cas9"
    tx_df <- count_df %>% filter(if_any(.cols = all_of(group_col),
                                        .fns = ~ .x == tx))
    tx_ns <- tx_df %>% pull(all_of(factor_col))
    missing_zero <- ns[!ns %in% tx_ns]
    if(length(missing_zero) > 0){
      zero_rows <- tx_df[1:length(missing_zero), ]
      zero_rows[, fill_col] <- 0
      zero_rows[, factor_col] <- missing_zero
      count_df <- rbind(count_df, zero_rows)
    }
  }
  count_df <- count_df %>% arrange(across(all_of(c(group_col, factor_col))))
  return(count_df)
}
