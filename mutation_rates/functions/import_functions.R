
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
  vcf_in$CHROM_ID <- factor(vcf_in$CHROM)
  chr_len <- length(levels(vcf_in$CHROM_ID))
  vcf_in$CHROM <- factor(vcf_in$CHROM_ID, labels = str_pad(1:chr_len, width = 2, pad = "0"))
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

HC_multiVcfParse_allVar <- function(vcf_file., all_Alts = F, gVCF = F, info_cols = NULL) {
  ## takes as input a VCF file output from GATK CombineGVCFs > GenotypeGVCFs
  # This file contains all samples from a founder group, with each sample's genotype
  # information occupying a single column
  # info_cols = c("QD", "FS", "MQ", "SOR")
  require(tidyverse)
  require(dplyr)
  require(reshape2)
  require(readr)
  # vcf_file. <- drive_const_file
  vcf_in <- read.table(vcf_file., stringsAsFactors = F)
  # vcf_in
  ## Extract column header from commented section of VCF
  # Pull first 100 lines from file as text
  hdr_full <- readLines(vcf_file., n = 120)
  
  contigs <- hdr_full[grep("^\\s*##contig=", hdr_full)]
  contig_lengths. = sapply(strsplit(contigs, ","), "[[", 2)
  chrom_lengths. <- gsub("[^0-9]", "", contig_lengths.) %>% as.numeric()
  # Line index of last "##" metadata header
  skip <- max(grep("^\\s*##", hdr_full))
  # Extract header line, tab-separated values to vector, remove "#" from CHROM, apply to dataframe
  hdr_raw <- hdr_full[skip+1]
  hdr_vctr <- unlist(strsplit(hdr_raw, split="\t"))
  hdr_vctr[1] <- substr(hdr_vctr[[1]], start = 2, stop=6)
  colnames(vcf_in) <- hdr_vctr
  # remove unused columns "ID" and "FILTER"
  vcf_parse <- subset(vcf_in, select = -c(ID, FILTER))
  rm(vcf_in)
  ## factor CHROM names, remove mitochondrial CHROM, make column of CHROM numbers
  vcf_parse$CHROM <- factor(vcf_parse$CHROM)
  vcf_parse <- vcf_parse %>% filter(CHROM != "NC_001224.1")
  vcf_parse$chrom_n <- as.numeric(vcf_parse$CHROM)
  vcf_parse$CHROM <- factor(str_pad(vcf_parse$chrom_n, width = 2, pad = "0"))
  vcf_parse$QUAL <- as.numeric(vcf_parse$QUAL)

  if(length(info_cols) > 0) {
    if(gVCF) {
      # INFO_col <- data.frame(INFO = c("End=100", "End=0", "BaseQ=5;MQ=10"))
      End_info <- subset(vcf_parse, select = c(CHROM, POS, INFO))
      colnames(End_info)[4] <- "temp"
      End_info <- cbind(End_info[, 1:2], colsplit(End_info$temp, pattern = "=", names = c("i", "End")))
      iEnd <- grep("END", End_info$i)
      End_info$End[-iEnd] <- NA
      End_info$End <- as.numeric(End_info$End)
      End_info$EndPOSi <- End_info$End
    } else{
      info_col <- subset(vcf_parse, select = c(POS, INFO))
      info_split <- colsplit(info_col$INFO, pattern = ";", names = paste0("V", 1:15))
      # extraInfo_rows <- apply(info_split, MARGIN = 1, function(x) length(grep("Rank", x)) != 0)
      # info_split[extraInfo_rows, ] <- info_split[extraInfo_rows, -c(4, 11, 13, 15)]
      na_cols <- apply(info_split, MARGIN = 2, function(x) sum(is.na(x), na.rm = T) > 1)
      info_split <- info_split[, !na_cols]
      info_names <- sapply(info_split[1,], function(x) str_split(x, pattern = "=", simplify = T)[1,1])
      colnames(info_split) <- info_names
      split_df <- data.frame(lapply(info_split, function(x) as.numeric(str_split(x, pattern = "=", simplify = T)[, 2])))
      info_df <- split_df[, info_cols]
    }
  }
  ## generate genomic position indicies
  vcf_parse$POSi <- vcf_parse$POS
  ## Chromosome lengths are from the S288C_R64_1.1 reference genome
  # chrom_lengths. <- c(230218, 813184, 316620, 1531933, 576874, 270161, 1090940, 562643, 439888, 745751,
  #                    666816, 1078177, 924431, 784333, 1091291, 948066)
  chrom_indcs <- cumsum(chrom_lengths.)
  # Calculate global positions based on chromosome positions and boundaries
  if(length(chrom_lengths.) >= 2) {
    for (chrom in 2:16) {
      # chrom <- 2
      #rm(chrom)
      chromNm <- levels(vcf_parse$CHROM)[chrom]
      idx <- vcf_parse$CHROM == chromNm
      vcf_parse[idx, "POSi"] <- vcf_parse$POS[idx] + chrom_indcs[chrom-1]
      if(gVCF) {
        End_info$EndPOSi[idx] <- End_info$End[idx] + chrom_indcs[chrom-1]
      }
    }
  }
  # vcf_parse$chrom_n <- as.numeric(vcf_parse$CHROM)
  # place POSi column next to POS
  vcf_parse <- cbind(subset(vcf_parse, select = c(CHROM, POS)), POSi = vcf_parse$POSi, 
                     subset(vcf_parse, select = -c(CHROM, POS, chrom_n, POSi)))
  ## transform FORMAT field into genotype, read depth, and genotype likelihood columns
  phsdColNames <- c("GT", "AlleleDP", "TotDP", "GQ", "PGT", "PID", "PhredLike", "PS")
  pos_cols <- subset(vcf_parse, select = c(CHROM, POS, POSi, REF, ALT, QUAL))
  if(length(info_cols) > 0) {
    pos_cols <- cbind(pos_cols, info_df)
  }
  # Loop through sample genotype columns to parse into dataframe columns
  vtmp <- subset(vcf_parse, select = -c(CHROM, POS, POSi, REF, ALT, QUAL, INFO, FORMAT))
  for (i in 1:ncol(vtmp)) {
    # i=1
    col_suffix <- colnames(vtmp)[c]
    FrmtCols <- colsplit(vtmp[,c], ":", phsdColNames)
    FrmtCols$ID <- substr(col_suffix, 1, 5)
    # in-phase genotypes are indicated by a "|", needs to be changed to "/" for consistency
    FrmtCols$GT <- gsub(pattern="|", replacement="/", FrmtCols$GT, fixed=T)
    
    ADCols <- colsplit(FrmtCols$AlleleDP, ",", c("Ref_DP", "Alt_DP"))
    ADCols$Ref_DP <- as.numeric(ADCols$Ref_DP)
    AltDPcols <- colsplit(ADCols$Alt_DP, ",", 
                          c("Alt1_DP", "Alt2_DP", "Alt3_DP", "Alt4_DP", "Alt5_DP", "Alt6_DP", "excess"))
    AltDPcols <- AltDPcols[, -ncol(AltDPcols)]
    AltDPcols <- as.data.frame(lapply(AltDPcols, as.numeric))
    if(gVCF) {
      AltDPcols[is.na(AltDPcols)] <- 0
    }
    # AltDPcols[is.na(AltDPcols)] <- 0
    fract_Ref_v <- ADCols$Ref_DP/(ADCols$Ref_DP + apply(AltDPcols, MARGIN = 1, sum, na.rm = T))
    Sum_DP_v <- ADCols$Ref_DP + apply(AltDPcols, MARGIN=1, sum, na.rm=T)
    
    sampleGTDP <- cbind(GT = FrmtCols$GT, 
                        Ref_DP = as.numeric(ADCols[, 1]), 
                        AltDPcols, 
                        Sum_DP = Sum_DP_v, 
                        fract_Ref = fract_Ref_v,
                        GQ = as.numeric(FrmtCols$GQ), 
                        ID = FrmtCols$ID)
    GTs <- unique(FrmtCols$GT)
    alt_GTs <- sort(GTs[!(GTs %in% c("./.", "0/0", "0/1", "1/1"))])
    GTlvls <- c("./.", "0/0", "0/1", "1/1", alt_GTs)
    sampleGTDP$GT <- factor(sampleGTDP$GT, levels=GTlvls)
    vcf_parse_sample <- cbind(pos_cols, sampleGTDP)
    print(paste0(FrmtCols$ID[1], " is done (", c, "/", ncol(vtmp), ")"))
    if (i == 1) {
      vcf_parsed <- vcf_parse_sample
    } else {
      vcf_parsed <- rbind(vcf_parsed, vcf_parse_sample)
    }
  }
  rm(vcf_parse)
  rm(vtmp)
  # Get rid of unused genotypes
  vcf_parsed$GT <- droplevels(vcf_parsed$GT)
  
  # Get rid of unused alternative allele columns
  # altIndx <- vcf_parsed %>% select(Alt1_DP:Alt6_DP) %>% apply(MARGIN = 2, sum, na.rm=T)
  AltDPcols <- subset(vcf_parsed, select = Alt1_DP:Alt6_DP)
  ##### Parse multiple alternative alleles
  AltCols <- colsplit(vcf_parsed$ALT, ",", c("ALT1", "ALT2", "ALT3", "ALT4", "ALT5", "ALT6", "excess"))
  AltCols <- subset(AltCols, select = -excess)
  AltCols <- as.data.frame(lapply(AltCols, as.character))
  AltCols <- replace(AltCols, AltCols == "", NA)
  AltCols <- AltCols[colSums(!is.na(AltCols)) > 0]
  if(gVCF) {
    vcf_parsed <- cbind(subset(vcf_parsed, select = CHROM:POSi), 
                        End = End_info$EndPOSi,
                        subset(vcf_parsed, select = REF:ID))
  }
  if (all_Alts == T) {
    vcf_parsed <- cbind(subset(vcf_parsed, select = CHROM:REF), 
                        AltCols,
                        subset(vcf_parsed, select = QUAL:Ref_DP), 
                        AltDPcols, 
                        subset(vcf_parsed, select = Sum_DP:ID))
  } else {
    vcf_parsed <- cbind(subset(vcf_parsed, select = CHROM:REF), 
                        ALT = AltCols[, 1],
                        subset(vcf_parsed, select = QUAL:Ref_DP), 
                        Alt_DP = AltDPcols[, 1], 
                        subset(vcf_parsed, select = Sum_DP:ID))
    vcf_parsed <- vcf_parsed %>% filter(GT %in% c("0/0", "0/1", "1/1"))
  }
  rm(AltCols)
  rm(AltDPcols)
  # Extract treatment, lineage, and replicate identifiers
  vcf_parsed$ID <- factor(vcf_parsed$ID)
  vcf_parsed$Tx <- factor(substr(vcf_parsed$ID, 1, 1))
  vcf_parsed$Line <- factor(substr(vcf_parsed$ID, 1, 3))
  vcf_parsed$Rep <- factor(substr(vcf_parsed$ID, 4, 5))
  
  return(vcf_parsed)
}

chainToDF <- function(chainFile) {
  chainTable <- read_lines(chainFile)
  chainTable <- gsub(pattern = "\"", replacement = "", x = chainTable)
  chainTable <- chainTable[!chainTable == ""]
  chrIdx <- which(chainTable %in% chrom_bound_BY$CHROM)
  endIdx <- length(chainTable)
  valueIdx <- data.frame(from = chrIdx + 1, to = c(chrIdx[2:16] - 1, endIdx))
  chrNames <- chainTable[chrIdx]
  
  chainDf <- data.frame(NULL)
  for (ch in 1:nrow(valueIdx)) {
    # ch = 2
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
    subDF$BY_POSi <- subDF$BY_POS + chrom_bound_BY$Start[ch] - 1
    subDF$RM_POSi <- subDF$RM_POS + chrom_bound_RM$Start[ch] - 1
    subDF$CHROM <- chrNames[ch]
    chainDf <- rbind(chainDf, subDF)
    chainDf$CHROM <- factor(chainDf$CHROM)
  }
  return(chainDf)
}

BYtoRMchainDf <- chainToDF(chainFile)

RMxBY_liftover <- function(vcf_df, chain_df = BYtoRMchainDf, 
                           lift_from = "RM", lift_to = "BY", 
                           chrom_from = chrom_bound_RM, chrom_to = chrom_bound_BY) {
  # vcf_df <- all_fltrd_RM %>% filter(ID == "F_A01")
  if(!(lift_to %in% c("BY", "RM"))) {
    print("lift_to must be 'BY' or 'RM'")
  } else {
    posi_col <- paste0(lift_from, "_POSi") 
    csDiff_col <- paste0("cs", lift_from, "diff")
  }
  chrom_diff <- c(0, chrom_from$End - chrom_to$End)
  vcf_df$POS_adj <- vcf_df$POS
  vcf_df$POSi_adj <- vcf_df$POSi
  for(ch in levels(chain_df$CHROM)) {
    # ch = levels(BYtoRMchainDf$CHROM)[6]
    ch_chain <- chain_df %>% filter(CHROM == ch)
    # iChrom_vcf <- vcf_df$CHROM == ch
    nBlocks <- nrow(ch_chain)
    for(b in 1:(nBlocks - 1)) {
      # b = 7
      b_i <- vcf_df$POSi >= ch_chain[b, posi_col] & vcf_df$POSi < ch_chain[(b + 1), posi_col]
      vcf_df$POS_adj[b_i] <- vcf_df$POS[b_i] + ch_chain[b, csDiff_col]
      vcf_df$POSi_adj[b_i] <- vcf_df$POSi[b_i] + ch_chain[b, csDiff_col] - chrom_diff[as.numeric(ch)]
    }
    e_i <- vcf_df$POSi >= ch_chain[nBlocks, posi_col] & vcf_df$POSi <= chrom_from$End[as.numeric(ch)]
    vcf_df$POS_adj[e_i] <- vcf_df$POS[e_i] + ch_chain[nBlocks, csDiff_col]
    vcf_df$POSi_adj[e_i] <- vcf_df$POSi[e_i] + ch_chain[nBlocks, csDiff_col] - chrom_diff[as.numeric(ch)]
    print(ch)
  }
  vcf_df$POS <- vcf_df$POS_adj
  vcf_df$POSi <- vcf_df$POSi_adj
  vcf_df <- vcf_df %>% dplyr::select(-c("POS_adj", "POSi_adj"))
  return(vcf_df)
}

import_depths <- function(depthDir, callSet = "BY") {
  # depthDir <- paste0(dataDir, "depth_metrics/")
  dir_list <- list.files(path = depthDir)
  depth_list <- dir_list[grep(callSet, dir_list)]
  depth_out <- data.frame(NULL)
  for (f in 1:length(depth_list)) {
    # f=1
    ID_path <- paste0(depthDir, depth_list[f])
    id <- substr(depth_list[f], 1, 5)
    tx <- substr(id, 1, 1)
    depth_table <- read_delim(ID_path, delim = "\t", 
                              skip = 8, comment = "#", progress = F)
    colnames(depth_table)[2] <- "count"
    nc_table <- depth_table %>% filter(coverage <= 6)
    total_nc <- sum(nc_table$count)
    total_sites <- sum(depth_table$count)
    stats_table <- read_delim(ID_path, delim = "\t", 
                              skip = 6, n_max = 1, progress = F) %>% 
      select(contains("MEDIAN"), contains("EXC"))
    ID_depth <- data.frame(ID = id, Tx = tx, callset = callSet, 
                           n_sites = total_sites, n_low_cover = total_nc, 
                           n_valid = total_sites - total_nc)
    ID_depth <- cbind(ID_depth, stats_table)
    depth_out <- rbind(depth_out, ID_depth)
  }
  return(depth_out)
}


# Unused ####################################################################################
HC_multiVcfParse_LOH <- function(vcf_file., refSequence, all_Alts = F, refVariants="", chrom_lengths.) {
  ## takes as input a VCF file output from GATK CombineGVCFs > GenotypeGVCFs
  # This file contains all samples from a founder group, with each sample's genotype
  # information occupying a single column
  require(tidyverse)
  require(dplyr)
  require(reshape2)
  require(readr)
  # vcf_file. <- paste0(allVarDir, "anc_pool_BYc_local.g.vcf")
  vcf_in <- read.table(vcf_file., stringsAsFactors = F)
  vcf_in
  ## Extract column header from commented section of VCF
  # Pull first 100 lines from file as text
  hdr_full <- readLines(vcf_file., n = 120)
  # Line index of last "##" metadata header
  skip <- max(grep("^\\s*##", hdr_full))
  # Extract header line, tab-separated values to vector, remove "#" from CHROM, apply to dataframe
  hdr_raw <- hdr_full[skip+1]
  hdr_vctr <- unlist(strsplit(hdr_raw, split="\t"))
  hdr_vctr[1] <- substr(hdr_vctr[[1]], start = 2, stop=6)
  colnames(vcf_in) <- hdr_vctr
  # remove unused columns "ID" and "FILTER"
  vcf_parse <- subset(vcf_in, select = -c(ID, FILTER))
  rm(vcf_in)
  ## factor CHROM names, remove mitochondrial CHROM, make column of CHROM numbers
  vcf_parse$CHROM <- factor(vcf_parse$CHROM)
  vcf_parse <- vcf_parse %>% filter(CHROM != "NC_001224.1")
  vcf_parse$chrom_n <- as.numeric(vcf_parse$CHROM)
  vcf_parse$CHROM <- factor(str_pad(vcf_parse$chrom_n, width = 2, pad = "0"))
  vcf_parse$QUAL <- as.numeric(vcf_parse$QUAL)
  ## generate genomic position indicies
  vcf_parse$POSi <- vcf_parse$POS
  ## Chromosome lengths are from the S288C_R64_1.1 reference genome
  # chrom_lengths. <- c(230218, 813184, 316620, 1531933, 576874, 270161, 1090940, 562643, 439888, 745751,
  #                    666816, 1078177, 924431, 784333, 1091291, 948066)
  chrom_indcs <- cumsum(chrom_lengths.)
  # Calculate global positions based on chromosome positions and boundaries
  for (chrom in 2:16) {
    #chrom<-2
    #rm(chrom)
    chromNm <- levels(vcf_parse$CHROM)[chrom]
    idx <- vcf_parse$CHROM == chromNm
    vcf_parse[idx, "POSi"] <- vcf_parse$POS[idx] + chrom_indcs[chrom - 1]
  }
  # vcf_parse$chrom_n <- as.numeric(vcf_parse$CHROM)
  # place POSi column next to POS
  vcf_parse <- cbind(subset(vcf_parse, select = c(CHROM, POS)), 
                     POSi = vcf_parse$POSi, 
                     subset(vcf_parse, select = -c(CHROM, POS, chrom_n, POSi)))
  
  ## transform FORMAT field into genotype, read depth, and genotype likelihood columns
  phsdColNames <- c("GT", "AlleleDP", "TotDP", "GQ", "PGT", "PID", "PhredLike", "PS")
  # Loop through sample genotype columns to parse into dataframe columns
  vtmp <- subset(vcf_parse, select = -c(CHROM, POS, POSi, REF, ALT, QUAL, INFO, FORMAT))
  for (i in 1:ncol(vtmp)) {
    # i=1
    col_suffix <- colnames(vtmp)[i]
    FrmtCols <- colsplit(vtmp[,i], ":", phsdColNames)
    FrmtCols$ID <- substr(col_suffix, 1, 5)
    # in-phase genotypes are indicated by a "|", needs to be changed to "/" for consistency
    FrmtCols$GT <- gsub(pattern="|", replacement="/", FrmtCols$GT, fixed=T)
    
    ADCols <- colsplit(FrmtCols$AlleleDP, ",", c("Ref_DP", "Alt_DP"))
    ADCols$Ref_DP <- as.numeric(ADCols$Ref_DP)
    AltDPcols <- colsplit(ADCols$Alt_DP, ",", 
                          c("Alt1_DP", "Alt2_DP", "Alt3_DP", "Alt4_DP", "Alt5_DP", "Alt6_DP"))
    AltDPcols <- as.data.frame(lapply(AltDPcols, as.numeric))
    AltDPcols[is.na(AltDPcols)] <- 0
    fract_Ref_v <- ADCols$Ref_DP/(ADCols$Ref_DP + AltDPcols[, 1])
    Sum_DP_v <- ADCols$Ref_DP + apply(AltDPcols, MARGIN = 1, sum, na.rm = T)
    
    sampleGTDP <- cbind(GT=FrmtCols$GT, Ref_DP=as.numeric(ADCols[, 1]), AltDPcols)
    sampleGTDP <- cbind(sampleGTDP, Sum_DP=Sum_DP_v, 
                        fract_Ref=fract_Ref_v,
                        GQ=as.numeric(FrmtCols$GQ), ID=FrmtCols$ID)
    GTlvls <- c("./.", "0/0", "0/1", "1/1", "0/2", "0/3", "0/4",
                "1/2", "1/3", "1/4", "2/2", "2/3", "2/4",
                "3/3", "3/4", "4/4")
    sampleGTDP$GT <- factor(sampleGTDP$GT, levels=GTlvls)
    vcf_parse_sample <- cbind(vcf_parse[, c(1:6)], sampleGTDP)
    if (i == 1) {
      vcf_parsed <- vcf_parse_sample
    } else {
      vcf_parsed <- rbind(vcf_parsed, vcf_parse_sample)
    }
  }
  # Get rid of unused genotypes
  vcf_parsed$GT <- droplevels(vcf_parsed$GT)
  
  # Get rid of unused alternative allele columns
  altIndx <- vcf_parsed %>% select(Alt1_DP:Alt6_DP) %>% apply(MARGIN = 2, sum, na.rm=T)
  AltDPcols <- vcf_parsed %>% select(Alt1_DP:Alt6_DP)
  ##### Parse multiple alternative alleles
  AltCols <- colsplit(vcf_parsed$ALT, ",", c("Alt1", "Alt2", "Alt3", "Alt4", "Alt5", "Alt6"))
  AltCols <- as.data.frame(lapply(AltCols, as.character))
  AltCols <- replace(AltCols, AltCols == "", NA)
  AltCols <- AltCols[colSums(!is.na(AltCols)) > 0]
  if (all_Alts == T) {
    vcf_parsed <- cbind(vcf_parsed[, c(1:4)], 
                        AltCols,
                        vcf_parsed[,6:8], 
                        AltDPcols[,altIndx > 0], 
                        vcf_parsed[,c(15:ncol(vcf_parsed))])
  } else {
    vcf_parsed <- cbind(vcf_parsed[,c(1:4)], 
                        ALT=AltCols[,1],
                        vcf_parsed[,6:8], 
                        Alt_DP=AltDPcols[,1], 
                        vcf_parsed[,c(15:ncol(vcf_parsed))])
    vcf_parsed <- vcf_parsed %>% filter(GT=="0/0" | GT=="0/1" | GT=="1/1")
  }
  # Extract treatment, lineage, and replicate identifiers
  vcf_parsed$ID <- factor(vcf_parsed$ID)
  vcf_parsed$Tx <- factor(substr(vcf_parsed$ID, 1, 1))
  vcf_parsed$Line <- factor(substr(vcf_parsed$ID, 1, 3))
  vcf_parsed$Rep <- factor(substr(vcf_parsed$ID, 4, 5))
  
  if (refVariants != "") {
    ## Merge vcf to reference variant file
    # refVariants <- RMxBYbcfVcf
    vcf_parsed <- merge(refVariants[,1:4], vcf_parsed, 
                        by=c("CHROM", "POS"), sort=T, 
                        all.x = T, all.y = T, suffixes = c(".ref", ".cln"))
  }
  return(vcf_parsed)
}

parseAllVCFs <- function(vcfPath) {
  # Produce list of VCF files in a folder, then applies "HC_multiVcfParse" 
  # to parse VCF files into a dataframe. Each founder group VCF is addended to the end
  # of the dataframe, such that each clone occupies a set of rows.
  require(DescTools)
  dir_list <- list.files(path=vcfPath)
  vcf_list <- dir_list[grep(".vcf", dir_list)]
  final_vcf <- data.frame(NULL)
  for (f in 1:length(vcf_list)) {
    # f=1
    line_path <- paste0(vcfPath, vcf_list[f])
    line_vcf <- HC_multiVcfParse_LOH(line_path)
    print(vcf_list[f])
    final_vcf <- rbind(final_vcf, line_vcf)
  }
  final_vcf$ID <- factor(final_vcf$ID)
  final_vcf$Tx <- factor(final_vcf$Tx)
  final_vcf$Line <- factor(final_vcf$Line)
  final_vcf$Rep <- factor(final_vcf$Rep)
  final_vcf$Tx_name <- final_vcf$Tx
  final_vcf$Tx_name <- Recode(final_vcf$Tx_name, "WT" = "N", "Cas9" = "H", "Drive" = "F")
  return(final_vcf)
}

import_GATKtable <- function(table_file, ccc = 1, refSeq = "BY") {
  ## takes as input a VCF file output from GATK CombineGVCFs > GenotypeGVCFs
  # This file contains all samples from a founder group, with each sample's genotype
  # information occupying a single column
  require(tidyverse)
  require(dplyr)
  require(reshape2)
  require(readr)
  
  if(refSeq == "BY") {
    chrom_bound. <- chrom_bound_BY
  } else if(refSeq == "RM") {
    chrom_bound. <- chrom_bound_RM
  } else{
    stop()
    print("reference must be RM or BY")
  }
  
  # table_file <- paste0(allVarDir, "/tables/A_BYaBYc_allVar.tsv")
  table_in <- read.table(table_file, stringsAsFactors = F, header = T)
  
  n_chrom <- length(unique(table_in$CHROM))
  table_in$CHROM <- factor(table_in$CHROM, labels = str_pad(1:n_chrom, 2, pad = "0"))
  table_in <- table_in %>% filter(CHROM %in% chrom_bound.$CHROM)
  
  site_cols <- grep(".", colnames(table_in), fixed = T, invert = T)
  sample_cols <- grep(".", colnames(table_in), fixed = T)
  table_sites <- table_in %>% select(colnames(table_in)[site_cols])
  table_samples <- table_in %>% select(colnames(table_in)[sample_cols])
  # rm(table_in)
  
  table_sites$POSi <- table_sites$POS
  # Calculate global positions based on chromosome positions and boundaries
  # Calculates POSi as POSi + cumsum(CHROM lengths) - 1
  for (chrom in chrom_bound.$CHROM) {
    #chrom<-2
    #rm(chrom)
    idx <- table_sites$CHROM == chrom
    chrom_diff <- chrom_bound.$Start[chrom_bound.$CHROM== chrom]
    table_sites[idx, "POSi"] <- table_sites$POS[idx] + chrom_diff - 1
  }
  
  # place POSi column next to POS
  table_sites <- cbind(subset(table_sites, select = c(CHROM, POS)), POSi = table_sites$POSi, 
                       subset(table_sites, select = -c(CHROM, POS, POSi)))
  
  table_samples <- cbind(POSi = table_sites$POSi, table_samples)
  gt_colnames <- lapply(str_split(colnames(table_samples)[2:11], fixed(".")), "[[", 2) %>% 
    unlist %>% unique()
  
  ## transform FORMAT field into genotype, read depth, and genotype likelihood columns
  
  for(c in 1:length(gt_colnames)) {
    gt_col <- table_samples %>% select(POSi, contains(gt_colnames[i])) %>%
      pivot_longer(cols = contains(gt_colnames[i]), names_to = c("ID", "discard"), 
                   names_sep = paste0(".", gt_colnames[i]), values_to = gt_colnames[i]) %>% 
      select(-discard)
    if(i == 1) {
      assemble_cols <- gt_col
    } else {
      assemble_cols <- merge(assemble_cols, gt_col, by = c("POSi", "ID"))
    }
    print(paste0(i, "/", length(gt_colnames), " done"))
  }
  rm(table_samples)
  assemble_cols <- assemble_cols %>% arrange(ID, POSi)
  
  assemble_cols$ID <- substr(assemble_cols$ID, 1, 5)
  # in-phase genotypes are indicated by a "|", needs to be changed to "/" for consistency
  assemble_cols$GT <- gsub(pattern="|", replacement="/", assemble_cols$GT, fixed=T)
  
  ADCols <- colsplit(assemble_cols$AD, ",", c("Ref_DP", "Alt_DP"))
  # ADCols$Ref_DP <- as.numeric(ADCols$Ref_DP)
  ADCols <- cbind(Ref_DP = ADCols$Ref_DP, colsplit(ADCols$Alt_DP, ",", 
                                                   "Alt_DP", "excess"))
  ADCols <- ADCols[, -ncol(ADCols)]
  i_drop_alt <- colSums(is.na(ADCols)) == nrow(ADCols)
  ADCols <- ADCols[, !i_drop_alt]
  ADCols <- as.data.frame(lapply(ADCols, as.numeric))
  if(gtable) {
    ADCols[is.na(ADCols)] <- 0
  }
  Sum_DP_v <- apply(ADCols, MARGIN=1, sum, na.rm=T)
  fract_Ref_v <- ADCols$Ref_DP/Sum_DP_v
  
  table_full <- cbind(table_sites,
                      GT = assemble_cols$GT,
                      ADCols, 
                      Sum_DP = Sum_DP_v, 
                      fract_Ref = fract_Ref_v,
                      GQ = as.numeric(assemble_cols$GQ), 
                      ID = assemble_cols$ID)
  
  rm(assemble_cols)
  rm(table_sites)
  rm(ADCols)
  # Parse Genotypes
  AltCols <- colsplit(table_full$ALT, ",", c("ALT", "excess"))
  # AltCols <- colsplit(table_full$ALT, ",", c(paste0("ALT", 1:n_Alts), "excess"))
  AltCols <- AltCols[, -ncol(AltCols)]
  i_drop_alt <- colSums(is.na(AltCols)) == nrow(AltCols)
  AltCols <- AltCols[, !i_drop_alt]
  GTcols <- colsplit(table_full$GT, "/", c("All_1", "All_2"))
  table_full$GT_crct <- "./."
  i_Ref_hom <- GTcols$All_1 == table_full$REF & GTcols$All_2 == table_full$REF
  i_Het_01 <- GTcols$All_1 == table_full$REF & GTcols$All_2 == AltCols$ALT1
  i_Alt_1 <- GTcols$All_1 == AltCols$ALT1 & GTcols$All_2 == AltCols$ALT1
  if(n_Alts > 1) {
    print("More alt alleles not yet supported")
  }
  table_full$GT_crct[i_Ref_hom] <- "0/0"
  table_full$GT_crct[i_Het_01] <- "0/1"
  table_full$GT_crct[i_Alt_1] <- "1/1"
  
  GTlvls <- unique(table_full$GT_crct) %>% sort()
  table_full$GT <- factor(table_full$GT_crct, levels = c(GTlvls[2:length(GTlvls)], "./."))
  
  table_full <- cbind(subset(table_full, select = CHROM:REF), 
                      ALT = AltCols$ALT,
                      subset(table_full, select = QUAL:ID))
  
  rm(AltCols)
  # 
  #   
  #   if(gtable) {
  #     table_parsed <- cbind(subset(table_parsed, select = CHROM:POSi), 
  #                         End = End_info$EndPOSi,
  #                         subset(table_parsed, select = REF:ID))
  #   }
  #   if (all_Alts == T) {
  #     table_parsed <- cbind(subset(table_parsed, select = CHROM:REF), 
  #                         AltCols,
  #                         subset(table_parsed, select = QUAL:Ref_DP), 
  #                         AltDPcols, 
  #                         subset(table_parsed, select = Sum_DP:ID))
  #   } else {
  #     table_parsed <- cbind(subset(table_parsed, select = CHROM:REF), 
  #                         ALT = AltCols[, 1],
  #                         subset(table_parsed, select = QUAL:Ref_DP), 
  #                         Alt_DP = AltDPcols[, 1], 
  #                         subset(table_parsed, select = Sum_DP:ID))
  #     table_parsed <- table_parsed %>% filter(GT %in% c("0/0", "0/1", "1/1"))
  #   }
  
  # Extract treatment, lineage, and replicate identifiers
  table_full$ID <- factor(table_full$ID)
  table_full$Tx <- factor(substr(table_full$ID, 1, 1))
  table_full$Line <- factor(substr(table_full$ID, 1, 3))
  table_full$Rep <- factor(substr(table_full$ID, 4, 5))
  
  return(table_full)
}


bcf_vcfParse <- function(vcf_path., clone_grp="all", refseq="BY") {
  # vcf_path. <- N_C_path
  if (clone_grp=="anc"){
    ptrn <- "00"
  } else if (clone_grp=="evo"){
    ptrn <- "[0-1][1-9]|10"
  } else if (clone_grp=="ND"){
    ptrn <- "N_"
  } else if (clone_grp=="HD"){
    ptrn <- "H_"
  } else if (clone_grp=="FD"){
    ptrn <- "F_"
  } else {
    ptrn <- ".vcf"
  }
  # vcf_path <- "~/SK_Lab/PhD_Projects/mutAccum/data/hcVariantCalls/RM_aligned/bcf/"
  # vcf_in <- read.table(paste0(vcf_path, "F_A01_RM_bcf.vcf"))
  path_list <- list.files(path=vcf_path., pattern=ptrn)
  vcf_all <- list()
  for (f in 1:length(path_list)) {
    # f=1
    # list of clone names for each dataframe
    line_nm <- substr(path_list[f],1,5)
    # read in vcf as a table and name columns
    vcf_in <- read.table(paste0(vcf_path., path_list[f]))
    VCFhdr <- c("CHROM", "POS", "ID", "REF", "ALT","QUAL", "FILTER", "INFO", "null", "FORMAT")
    colnames(vcf_in) <- VCFhdr
    # View(vcf_in)
    vcf_in$CHROM <- factor(vcf_in$CHROM)
    chrom_lengths <- c(230218, 813184, 316620, 1531933, 576874, 270161, 1090940, 562643, 439888, 745751,
                       666816, 1078177, 924431, 784333, 1091291, 948066)
    chrom_indcs <- cumsum(chrom_lengths)
    vcf_in$POSi <- vcf_in$POS
    for (chrom in 2:16) {
      #chrom<-2
      #rm(chrom)
      chromNm <- levels(vcf_in$CHROM)[chrom]
      idx <- vcf_in$CHROM==chromNm
      vcf_in[idx,"POSi"] <- vcf_in$POS[idx] + chrom_indcs[chrom-1]
    }
    # Parse the INFO filed to get depth and quality data
    vcf_in$Tot_DP <- -1
    vcf_in$Ref_DP <- -1
    vcf_in$Alt_DP <- -1
    vcf_in$Sum_DP <- -1
    vcf_in$MQ <- -1
    vcf_in$fract_BY <- -0.5
    
    str_INFO_split <- mapply(str_split, vcf_in$INFO, MoreArgs = list(pattern=";"), SIMPLIFY = T)
    
    vcf_in$Tot_DP <- sapply(str_INFO_split, grep, pattern="DP=", value=T) %>%
      sapply(., str_split, pattern="=", simplify=F) %>% unlist() %>% matrix(., ncol=2, byrow=T) %>% 
      subset(., select = 2) %>% as.numeric()
    
    DP4_vals <- sapply(str_INFO_split, grep, pattern="DP4=", value=T) %>% 
      sapply(., str_split, pattern="=", simplify=F) %>% unlist() %>% 
      matrix(., ncol=2, byrow=T) %>% subset(., select = 2) %>% 
      lapply(., str_split, pattern=",", simplify=F) %>% 
      unlist() %>% as.numeric() %>% matrix(., ncol=4, byrow=T)
    
    vcf_in$Ref_DP <- DP4_vals[,1] + DP4_vals[,2]
    vcf_in$Alt_DP <- DP4_vals[,3] + DP4_vals[,4]
    
    vcf_in$Sum_DP <- vcf_in$Ref_DP + vcf_in$Alt_DP
    
    vcf_in$MQ <- sapply(str_INFO_split, grep, pattern="MQ=", value=T) %>%
      sapply(., str_split, pattern="=", simplify=F) %>% unlist() %>% matrix(., ncol=2, byrow=T) %>% 
      subset(., select = 2) %>% as.numeric()
    
    vcf_in$fract_BY <- vcf_in$Ref_DP/(vcf_in$Ref_DP+vcf_in$Alt_DP)
    
    
    # parse the FORMAT field into genotypes and likelihoods
    GTcolNames <- c("GT", "LogL")
    GTcols <- colsplit(vcf_in$FORMAT, ":", GTcolNames)
    vcf_in$GT <- GTcols[ ,1]
    vcf_in$GT <- as.factor(vcf_in$GT)
    #*********#
    # If likelihoods are needed, can finish the rest of this code
    # lkCols <- colsplit(GTcols[,2], ",", names=c("ref","het","alt"))
    # lkCols[is.na(lkCols)] <- -1
    # GTlklh <- lkCols
    # for (l in 1:nrow(lkCols)) {
    #   GTlklh[l,] <- sort(lkCols[l,])
    # }
    # GTlklh <- t(sapply(t(lkCols), sort)) %>% matrix(., ncol=3, byrow=T) %>% subset(., select=2)
    # vcf_in$lklhd <- as.numeric(GTlklh)
    #********#
    vcf_in <- vcf_in %>% drop_na(GT)
    vcf_parsed <- vcf_in[,c(1,2,4,5,6,11:18)]
    # View(vcf_parsed)
    if (refseq=="RM") {
      vcf_parsed$REF <- vcf_in$Alt
      vcf_parsed$ALT <- vcf_in$Ref
      vcf_parsed$REF_DP <- vcf_in$Alt_DP
      vcf_parsed$ALT_DP <- vcf_in$Ref_DP
      levels(vcf_parsed$GT) <- c("NC", "RM", "het", "BY", "mut")
    } else if (refseq=="BY"){
      nGTs <- length(levels(vcf_parsed$GT))
      # levels(vcf_parsed$GT) <- c("NC", "BY", "het", "RM", paste0("mut", 1:(nGTs-4)))
    }
    else {
      print("reference not recognized")
    }
    vcf_parsed[,c(2,5,7,8,9)] <- sapply(vcf_parsed[,c(2,5,7,8,9)], as.numeric)
    
    vcf_all[[line_nm]] <- vcf_parsed
  }
  return(vcf_all)
}

bcf_multiVcfParse <- function(vcf_file., refSequence, all_Alts=F, refVariants="") {
  
  # read in vcf as a table and name columns
  # vcf_file. <- "~/SK_Lab/PhD_Projects/mutAccum/data/bcfCall/N_C_BYm.b.vcf"
  vcf_in <- read.table(vcf_file., stringsAsFactors = F)
  ## Extract column header from commented section of VCF
  # Pull first 100 lines from file as text
  hdr_full <- readLines(vcf_file., n = 100)
  # Line index of last "##" metadata header
  skip <- max(grep("^\\s*##", hdr_full))
  # Extract header line, tsv to vector, remove "#" from CHROM, apply to dataframe
  hdr_raw <- hdr_full[skip+1]
  hdr_vctr <- unlist(strsplit(hdr_raw, split="\t"))
  hdr_vctr[1] <- substr(hdr_vctr[[1]], start = 2, stop=6)
  colnames(vcf_in) <- hdr_vctr
  # View(vcf_in)
  vcf_in <- vcf_in %>% filter(CHROM != "NC_001224.1")
  vcf_in$CHROM <- factor(vcf_in$CHROM)
  
  chrom_lengths <- c(230218, 813184, 316620, 1531933, 576874, 270161, 1090940, 562643, 439888, 745751,
                     666816, 1078177, 924431, 784333, 1091291, 948066)
  chrom_indcs <- cumsum(chrom_lengths)
  vcf_in$POSi <- vcf_in$POS
  for (chrom in 2:16) {
    #chrom<-2
    #rm(chrom)
    chromNm <- levels(vcf_in$CHROM)[chrom]
    idx <- vcf_in$CHROM==chromNm
    vcf_in[idx,"POSi"] <- vcf_in$POS[idx] + chrom_indcs[chrom-1]
  }
  
  vcf_parse <- vcf_in[,c(1,2,ncol(vcf_in), 4,5,6,8,10:(ncol(vcf_in)-1))]
  vcf_parse$CHROM <- factor(vcf_parse$CHROM)
  vcf_parse$QUAL <- as.numeric(vcf_parse$QUAL)
  
  phsdColNames <- c("GT", "PhredLike")
  # Loop through sample genotype columns to parse into dataframe columns
  vtmp <- vcf_parse[,8:ncol(vcf_parse)]
  for (i in 1:ncol(vtmp)) {
    # i=1
    col_suffix <- colnames(vtmp)[c]
    FrmtCols <- colsplit(vtmp[,c], ":", phsdColNames)
    FrmtCols$ID <- substr(col_suffix, 1, 5)
    likeNames <- c("BYlike", "hetLike", "RMlike")
    likeCols <- colsplit(FrmtCols$PhredLike, ",", likeNames)
    likeCols_n <- as.data.frame(apply(likeCols, MARGIN = 2, as.numeric))
    GTlikeCols <- cbind(FrmtCols[,c(3,1)], likeCols_n)
    GTlikeCols$GTlike <- GTlikeCols$hetLike
    
    BYhomLikes <- GTlikeCols[GTlikeCols$GT=="0/0",c(4:5)]
    BYhomLike <- apply(BYhomLikes, MARGIN = 1, min)
    GTlikeCols$GTlike[GTlikeCols$GT=="0/0"] <- BYhomLike
    
    RMhomLikes <- GTlikeCols[GTlikeCols$GT=="1/1",c(3:4)]
    RMhomLike <- apply(RMhomLikes, MARGIN = 1, min)
    GTlikeCols$GTlike[GTlikeCols$GT=="1/1"] <- RMhomLike
    
    hetLikes <- GTlikeCols[GTlikeCols$GT=="0/1",c(3,5)]
    hetLike <- apply(hetLikes, MARGIN = 1, min)
    GTlikeCols$GTlike[GTlikeCols$GT=="0/1"] <- hetLike
    
    if (i == 1) {
      vcf_parsed <- cbind(vcf_parse[,c(1:6)], sampleGTDP)
    } else {
      vcf_parsed <- rbind(vcf_parsed, sampleGTDP)
    }
  }
  # Parse the INFO filed to get depth and quality data
  vcf_in$TOT_DP <- -1
  vcf_in$REF_DP <- -1
  vcf_in$ALT_DP <- -1
  vcf_in$SUM_DP <- -1
  vcf_in$MQ <- -1
  vcf_in$fract_BY <- -0.5
  
  str_INFO_split <- mapply(str_split, vcf_in$INFO, MoreArgs = list(pattern=";"), SIMPLIFY = T)
  
  vcf_in$TOT_DP <- sapply(str_INFO_split, grep, pattern="DP=", value=T) %>%
    sapply(., str_split, pattern="=", simplify=F) %>% unlist() %>% matrix(., ncol=2, byrow=T) %>% 
    subset(., select = 2) %>% as.numeric()
  
  DP4_vals <- sapply(str_INFO_split, grep, pattern="DP4=", value=T) %>% 
    sapply(., str_split, pattern="=", simplify=F) %>% unlist() %>% 
    matrix(., ncol=2, byrow=T) %>% subset(., select = 2) %>% 
    lapply(., str_split, pattern=",", simplify=F) %>% 
    unlist() %>% as.numeric() %>% matrix(., ncol=4, byrow=T)
  
  vcf_in$REF_DP <- DP4_vals[,1] + DP4_vals[,2]
  vcf_in$ALT_DP <- DP4_vals[,3] + DP4_vals[,4]
  
  vcf_in$SUM_DP <- vcf_in$REF_DP + vcf_in$ALT_DP
  
  vcf_in$MQ <- sapply(str_INFO_split, grep, pattern="MQ=", value=T) %>%
    sapply(., str_split, pattern="=", simplify=F) %>% unlist() %>% matrix(., ncol=2, byrow=T) %>% 
    subset(., select = 2) %>% as.numeric()
  
  vcf_in$fract_BY <- vcf_in$REF_DP/(vcf_in$REF_DP+vcf_in$ALT_DP)
  
  
  # parse the FORMAT field into genotypes and likelihoods
  GTcolNames <- c("GT", "LogL")
  GTcols <- colsplit(vcf_in$FORMAT, ":", GTcolNames)
  vcf_in$GT <- GTcols[ ,1]
  vcf_in$GT <- as.factor(vcf_in$GT)
  #*********#
  # If likelihoods are needed, can finish the rest of this code
  # lkCols <- colsplit(GTcols[,2], ",", names=c("ref","het","alt"))
  # lkCols[is.na(lkCols)] <- -1
  # GTlklh <- lkCols
  # for (l in 1:nrow(lkCols)) {
  #   GTlklh[l,] <- sort(lkCols[l,])
  # }
  # GTlklh <- t(sapply(t(lkCols), sort)) %>% matrix(., ncol=3, byrow=T) %>% subset(., select=2)
  # vcf_in$lklhd <- as.numeric(GTlklh)
  #********#
  vcf_in <- vcf_in %>% drop_na(GT)
  vcf_parsed <- vcf_in[,c(1,2,4,5,6,11:17)]
  # View(vcf_parsed)
  if (refseq=="RM") {
    vcf_parsed$REF <- vcf_in$ALT
    vcf_parsed$ALT <- vcf_in$REF
    vcf_parsed$REF_DP <- vcf_in$ALT_DP
    vcf_parsed$ALT_DP <- vcf_in$REF_DP
    levels(vcf_parsed$GT) <- c("NC", "RM", "het", "BY", "mut")
  } else if (refseq=="BY"){
    levels(vcf_parsed$GT) <- c("BY", "het", "RM", "mut")
  }
  else {
    print("reference not recognized")
  }
  vcf_parsed[,c(2,5,7,8,9)] <- sapply(vcf_parsed[,c(2,5,7,8,9)], as.numeric)
  
  vcf_all[[line_nm]] <- vcf_parse
  return(vcf_all)
}

HC_vcfParse <- function(vcf_path., clone_grp, refseq) {
  # vcf_path. <- NC00_HC_filenm
  if (clone_grp=="all"){
    ptrn <- ".vcf"
  } else if (clone_grp=="anc"){
    ptrn <- "00"
  } else if (clone_grp=="evo"){
    ptrn <- "[0-1][1-9]|10"
  } else if (clone_grp=="ND"){
    ptrn <- "N"
  }  else if (clone_grp=="HD"){
    ptrn <- "H"
  } else if (clone_grp=="FD"){
    ptrn <- "F"
  }
  path_list <- list.files(path=vcf_path., pattern=ptrn)
  vcf_all <- list()
  for (f in 1:length(path_list)) {
    # f=1
    line_nm <- substr(path_list[f],1,5)
    vcf_in <- read.table(paste0(vcf_path., path_list[f]))
    if(refseq=="BY"){
      VCFhdrBY <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "GT")
      colnames(vcf_in) <- VCFhdrBY
      newColNames <- c("GT", "RefAltDP", "TotDP", "GQ", "PhredLike")
      newCols <- colsplit(vcf_in[,10], ":", newColNames)
      
      newColsBY <- colsplit(newCols[,2], ",", c("Ref_DP", "Alt_DP"))
      allNewCols <- cbind(GT=newCols[,1], newColsBY, newCols[,3:4])
      allNewCols$GT <- factor(allNewCols$GT)
      allNewCols$Alt_DP <- as.numeric(allNewCols$Alt_DP)
      # levels(allNewCols$GT) <- c("BY", "het", "RM", "mut")
      vcf_parsed <- cbind(vcf_in[,c(1,2,4:6)], allNewCols)
    } else if(refseq=="RM"){
      VCFhdrRM <- c("CHROM", "POS", "ID", "ALT", "REF", "QUAL", "FILTER", "INFO", "FORMAT", "GT")
      colnames(vcf_in) <- VCFhdrRM
      # View(vcf_in)
      
      newColNames <- c("GT", "RefAltDP", "TotDP", "GQ", "PhredLike")
      newCols <- colsplit(vcf_in[,10], ":", newColNames)
      newColsRM <- colsplit(newCols[,2], ",", c("Alt_DP", "Ref_DP"))
      allNewCols <- cbind(GT=newCols[,1], newColsRM[,c(2,1)], newCols[,3:4])
      allNewCols$GT <- factor(allNewCols$GT, levels=c("0/0", "0/1", "1/1", "1/2"))
      levels(allNewCols$GT) <- c("RM", "het", "BY", "mut")
      vcf_parsed <- cbind(vcf_in[,c(1,2,5,4,6)], allNewCols)
    } else {print("reference not recognized")}
    
    chrom_lengths <- c(230218, 813184, 316620, 1531933, 576874, 270161, 1090940, 562643, 439888, 745751,
                       666816, 1078177, 924431, 784333, 1091291, 948066)
    chrom_indcs <- cumsum(chrom_lengths)
    
    # vcf_parse <- vcf_in[,-c(3,7)]
    vcf_parsed$CHROM <- factor(vcf_parsed$CHROM)
    vcf_parsed <- vcf_parsed %>% filter(CHROM != "NC_001224.1")
    vcf_parsed$QUAL <- as.numeric(vcf_parsed$QUAL)
    
    vcf_parsed$POSi <- vcf_parsed$POS
    for (chrom in 2:16) {
      #chrom<-2
      #rm(chrom)
      chromNm <- levels(vcf_parsed$CHROM)[chrom]
      idx <- vcf_parsed$CHROM==chromNm
      vcf_parsed[idx,"POSi"] <- vcf_parsed$POS[idx] + chrom_indcs[chrom-1]
    }
    str(vcf_parsed)
    vcf_all[[line_nm]] <- vcf_parsed
  }
  return(vcf_all)
}

