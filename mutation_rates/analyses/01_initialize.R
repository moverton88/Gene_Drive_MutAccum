
####################################################################################
# R script for the analysis of mutAccum GATK VCFs according to lineage
library(boot)
library(patchwork)
library(htmlTable)
library(kableExtra)
library(magick)
library(reshape2)
library(IRanges)
library(seqinr)
library(plyr)
library(pammtools)
library(slider)
library(cowplot)
# library(ggplot2)
library(RColorBrewer)
library(ggrepel)
library(coin)
library(ggpubr)
library(DescTools)
library(lsr)
library(rcompanion)
library(permute)
library(gridExtra)
library(data.table)
library(bedr, verbose = T)
library(kSamples)
library(uniftest)
library(effsize)
library(foreach)
library(twosamples)
library(doParallel)
library(scales)
library(gt)
library(tidyverse)

####################################################################################
# Set directory variables
setwd("/Users/Mastermind/SK_Lab/PhD_Projects/geneDrive/")

# Set directory paths
codeDir <- "./code/mutation_rates/"
dataDir <- "./data/"
dataIntDir <- "./data/int/"
posDir <- "./refseq/POS_files/"
allVarDir <- "./data/allVarVcfs/"
outDir <- "./results/"
outIntDir <- "./results/int/"

# Filenames for loading objects
chainFile <- "./refseq/RM_ref/variant_POS/RMxBY_chain.txt"
anc_file <- paste0(allVarDir, "anc_pool_BYc_local.vcf")
BY_file <- paste0(allVarDir, "A_BYaBYc_SNPs_filtered.vcf")
RM_file <- paste0(allVarDir, "A_RMaRMc_SNPs_filtered.vcf")
BY_indel_file <- paste0(allVarDir, "indels/", "A_BYaBYc_indel.vcf")
RM_indel_file <- paste0(allVarDir, "indels/", "A_RMaRMc_indel.vcf")
RMxBY_file <- file.path(posDir, "RMxBY_ref_noMit.vcf")
RMxBY_comp_file <- file.path(posDir, "RMxBY_comp_ref_bcf.vcf")
repeats_file <- paste0(posDir, "BY_repetitiveElements.bed")
colony_fn <- paste0(dataDir, "colony_size/2019.06.07_colony_sizes.tsv")
D_drive_const_file <- paste0(allVarDir, "Drive/F_Cas9_default_allVar.vcf")
C_drive_const_file <- paste0(allVarDir, "Drive/H_Cas9_default_allVar.vcf")
W_drive_const_file <- paste0(allVarDir, "Drive/N_Cas9_default_allVar.vcf")
# Sui_GC_file <- paste0(dataDir, "Sui_GC_lengths.txt")

# Filenames for saving objects
depth_filename <- paste0(dataIntDir, "seq_depths.RData")
BY_raw_filename <- paste0(dataIntDir, "BY_raw_2022_03.RData")
RM_raw_filename <- paste0(dataIntDir, "RM_raw_2022_03.RData")
SNPs_merge_raw_filename <- paste0(dataIntDir, "SNPs_merge_raw_2022_03.RData")
SNPs_merge_filename <- paste0(dataIntDir, "SNPs_merge_final_2022_03.RData")
LOH_SNPs_file <- paste0(dataIntDir, "LOH_SNPs_2022_03.RData")
denovo_SNPs_file <- paste0(dataIntDir, "denovo_SNPs_2022_03.RData")
bounds_filename <- paste0(dataIntDir, "all_LOHbounds_NS_2022_04.RData")
countsEC_filename <- paste0(dataIntDir, "all_LOHcounts_NS_2022_04.RData")
counts_filename <- paste0(dataIntDir, "all_LOHcounts_2022_04.RData")

triploid_merge_raw_filename <- paste0(dataIntDir, "triploid_SNPs_merge.RData")
triploid_merge_filename <- paste0(dataIntDir, "triploid_SNPs_final.RData")
triploid_LOH_SNPs_file <- paste0(dataIntDir, "triploid_LOH_SNPs.RData")
triploid_denovo_SNPs_file <- paste0(dataIntDir, "triploid_denovo_SNPs.RData")

# Define options and asthetics
options(scipen = 3, max.print = 100)

theme_set(theme_minimal() + 
            theme(plot.background = element_rect(fill = "white", colour = NA)))

lineagePal <-c(brewer.pal(8,"Dark2"), brewer.pal(4,"Set1"), brewer.pal(8,"Set2"))
txPal <- c("#444444", "#BAD600", "#40B0A6")
allelePal <- c("orange2", "grey20", "blue3")

GT_levels <- c("0/0", "0/1", "1/1")
Tx_levels <- c("N", "H", "F")
Tx_ID_levels <- c("W", "C", "D")
Tx_name_levels <- c("WT", "Cas9", "Drive")
Tx_combo <- data.frame(Tx_1 = c("WT", "WT", "Cas9"), 
                       Tx_2 = c("Cas9", "Drive", "Drive"))
Tx_ID_combo <- data.frame(Tx_1 = c("W", "W", "C"), 
                          Tx_2 = c("C", "D", "D"))
# End-point clones tested for drive activity
DA_clones <- c("F_A09", "F_C02", "F_D01", "F_F03", "F_F07", "F_G10")

# Chormosome indexing objects ------
roman_chr <- factor(as.character(as.roman(1:16)), levels = as.character(as.roman(1:16)))
chrom_IDs <- data.frame(num = 1:16, CHROM = str_pad(1:16, 2, pad = "0"), rom_CHROM = roman_chr,
                        ID = c("NC_001133.9", "NC_001134.8", "NC_001135.5", "NC_001136.10",
                               "NC_001137.3", "NC_001138.5", "NC_001139.9", "NC_001140.6",
                               "NC_001141.2", "NC_001142.9", "NC_001143.9", "NC_001144.5",
                               "NC_001145.3", "NC_001146.8", "NC_001147.6", "NC_001148.4"))

chrom_lengths_BY <- c(230218, 813184, 316620, 1531933, 576874, 270161, 1090940,
                      562643, 439888, 745751, 666816, 1078177, 924431, 784333, 
                      1091291, 948066)


chrom_lengths_BY_df <- cbind(chrom_IDs, chrom_length = chrom_lengths_BY)
chrom_lengths_BY_df$CHROM <- factor(chrom_lengths_BY_df$CHROM)

chrom_lengths_RM <- c(230348, 812918, 316506, 1532315, 576771, 270224, 1091128,
                      563026, 440151, 745731, 667107, 1078742, 924979, 784407,
                      1091989, 948187)

chrom_indcs <- cumsum(chrom_lengths_BY)

chrom_coord <- function(lngths) {
  chr_start <- c(1)
  for (l in 1:(length(lngths)-1)) {
    l_sum <- sum(lngths[1:l])+1
    chr_start <- c(chr_start, l_sum)
  }
  chr_end <- c()
  for (l in 1:length(lngths)) {
    l_sum <- sum(lngths[1:l])
    chr_end <- c(chr_end, l_sum)
  }
  chr_bound <- data.frame(Start=as.numeric(chr_start), 
                          End=as.numeric(chr_end), 
                          CHROM=factor(str_pad(1:length(lngths), width = 2, pad = "0")), 
                          stringsAsFactors = F)
  return(chr_bound)
}

chrom_bound_BY <- chrom_coord(chrom_lengths_BY)
chrom_bound_BY$length <- chrom_lengths_BY_df$chrom_length
chrom_bound_BY$rom_CHROM <- roman_chr[as.numeric(chrom_bound_BY$CHROM)]
chrom_bound_RM <- chrom_coord(chrom_lengths_RM)


centrom_POS <- c(151515, 238255, 114435, 449760,
                 152040, 148560, 496985, 105650,
                 355680, 436375, 440195, 150895,
                 268085, 628810, 326650, 556010)
centrom_POSi <- centrom_POS + chrom_bound_BY$Start - 1


centrom_df <- data.frame(CHROM = factor(str_pad(1:16, width = 2, pad = "0")), 
                         POS = centrom_POS, POSi = centrom_POSi)
centrom_df$rom_CHROM <- roman_chr[as.numeric(centrom_df$CHROM)]

chrom_arms <- data.frame(CHROM = centrom_df$CHROM, 
                         arm_1_cent = -centrom_df$POS, 
                         arm_2_cent = chrom_lengths_BY - centrom_df$POS)

g_length <- sum(chrom_lengths_BY)

# Problematic regions

repeats_bed <- read_delim(repeats_file, delim = "\t", col_names = c("CHROM", "Start_POS", "End_POS", "Name"))
repeats_bed <- repeats_bed %>% filter(CHROM <= 16)
repeats_bed$CHROM <- str_pad(repeats_bed$CHROM, 2, pad = "0")
repeats_bed$Start_POS <- repeats_bed$Start_POS + 1
repeats_bed$Start_POSi <- repeats_bed %>% ConvertPosIndicies(pos_col = "Start_POS", index_out = "POSi")
repeats_bed$End_POSi <- repeats_bed %>% ConvertPosIndicies(pos_col = "End_POS", index_out = "POSi")


rRNA_POSi <- data.frame(ID = "rRNA", POS = c(447000, 490000), 
                        POSi = chrom_bound_BY$Start[12] + c(447000, 490000))
Ty_POSi <- data.frame(ID = "Ty", POS = c(197400, 197600), 
                      POSi = chrom_bound_BY$Start[10] + c(197400, 197600))

# Create liftover file for RM vcf import
###############################################################################
BYtoRMchainDf <- chainToDF(chainFile)

# RMxBY variant site positions
###############################################################################

# Table created by aligning RM and BY references and outputting variant sites
# using bcftools consensus
RMxBY_vcf <- RMxBY_vcfParse(RMxBY_file)
iRMxBY_SNPs <- nchar(RMxBY_vcf$BY) == 1 & nchar(RMxBY_vcf$RM) == 1 
RMxBY_SNPs_POSi <- RMxBY_vcf$POSi[iRMxBY_SNPs]

iRMxBY_indels <- nchar(RMxBY_vcf$BY) > 1 | nchar(RMxBY_vcf$RM) > 1 
RMxBY_indels_POSi <- RMxBY_vcf$POSi[iRMxBY_indels]

RMxBY_comp_vcf <- RMxBY_vcfParse(RMxBY_comp_file)
iRMxBY_comp_SNPs <- nchar(RMxBY_comp_vcf$BY) == 1 & nchar(RMxBY_comp_vcf$RM) == 1 
RMxBY_comp_SNPs_POSi <- RMxBY_comp_vcf$POSi[iRMxBY_comp_SNPs]

# These sites do not show a consistent genotype in the ancestors
# potentially due to mapping difficulty or mutations
prob_sites <- c(187555, 187558, 3675443, 3675444, 3675448)

# Number of generations
tx_popGens <- data.frame(Tx_name = c("WT", "Cas9", "Drive"), 
                         pop_n = c(445200, 284545, 250027), 
                         gens = c( 826, 797, 788), stringsAsFactors = T)

n_gens <- mean(tx_popGens$gens)


BY_ref_fn <- "./refseq/BY_ref/S288C_R64_refseq.fasta"
BY_seq <- read.fasta(BY_ref_fn, as.string = T)
BY_seq_v <- read.fasta(BY_ref_fn, as.string = F)

BY4742_ref_fn <- "./refseq/BY_ref/BY4742_refseq.fasta"
BY4742_seq <- read.fasta(BY4742_ref_fn, as.string = T)
BY4742_seq_v <- read.fasta(BY4742_ref_fn, as.string = F)

RM_ref_fn <- "./refseq/RM_ref/RM_refseq_UCSD_2020_v4.fna"
RM_seq <- read.fasta(RM_ref_fn, as.string = T)
RM_seq_v <- read.fasta(RM_ref_fn, as.string = F)

gRNA_seq <- "acttgaagattctttagtgt"
gRNA_seq_split <- unlist(strsplit(gRNA_seq, split = ""))
gRNA_seq_comp_split <- rev(comp(gRNA_seq_split))
gRNA_seq_comp <- paste(gRNA_seq_comp_split, collapse = '')
target_POS <- list("start" = 565935, "end" = 565954)
