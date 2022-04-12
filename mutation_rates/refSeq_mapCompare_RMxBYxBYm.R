


old_theme <- theme_get()
theme_set(theme_minimal())

setwd("/Users/Mastermind/SK_Lab/PhD_Projects/mutAccum/code/R_MA_seq")

RM_vcf="~/SK_Lab/PhD_Projects/mutAccum/data/hcVariantCalls/RM_aligned/RMxBY_sites/anc_RM_RMxBY.vcf"
BY_vcf="~/SK_Lab/PhD_Projects/mutAccum/data/hcVariantCalls/BY_aligned/RMxBY_sites/anc_BY_RMxBY.vcf"
BYm_vcf="~/SK_Lab/PhD_Projects/mutAccum/data/hcVariantCalls/BYm_aligned/anc_BYm_RMxBY.vcf"

anc12_RM <- HC_multiVcfParse(RM_vcf, refSequence="BY")
anc12_BY <- HC_multiVcfParse(BY_vcf, refSequence="BY")
anc12_BYm <- HC_multiVcfParse(BYm_vcf, refSequence="BY")


anc12_RM_hiQual <- anc12_RM %>% filter(GenQual>90) %>% filter(GT %in% c("0/0", "0/1", "1/1"))
anc12_RM_hiQual$algn <- "RM"
nrow(anc12_RM_hiQual[anc12_RM_hiQual$ID == "F_A00",])

anc12_BY_hiQual <- anc12_BY %>% filter(GenQual>90) %>% filter(GT %in% c("0/0", "0/1", "1/1"))
anc12_BY_hiQual$algn <- "BY"
nrow(anc12_BY_hiQual[anc12_BY_hiQual$ID == "F_A00",])
nrow(anc12_BY_hiQual[anc12_BY_hiQual$ID == "N_E00",])

anc12_BYm_hiQual <- anc12_BYm %>% filter(GenQual>90) %>% filter(GT %in% c("0/0", "0/1", "1/1"))
anc12_BYm_hiQual$algn <- "BYm"
nrow(anc12_BYm_hiQual[anc12_BYm_hiQual$ID == "F_A00",])

anc_all_hiQual <- rbind(anc12_RM_hiQual, anc12_BY_hiQual, anc12_BYm_hiQual)

df1 <- anc12_BYm %>% filter(GenQual>90) %>%
  dplyr::count(ID, GT) %>%
  dplyr::group_by(ID) %>% #change to `group_by(Genotypes) %>%` for alternative approach
  mutate(prop = n / sum(n), count=n)
df1$algn <- "BYm"

df2 <- anc12_BY %>% filter(GenQual>90) %>%
  dplyr::count(ID, GT) %>%
  dplyr::group_by(ID) %>% #change to `group_by(Genotypes) %>%` for alternative approach
  mutate(prop = n / sum(n), count=n)
df2$algn <- "BY"

df3 <- anc12_RM %>% filter(GenQual>90) %>%
  dplyr::count(ID, GT) %>%
  dplyr::group_by(ID) %>% #change to `group_by(Genotypes) %>%` for alternative approach
  mutate(prop = n / sum(n), count=n)
df3$algn <- "RM"

df4 <- rbind(df1, df2, df3)

# Only genotypes 0/0 0/1 1/1
df4$GT <- droplevels(df4$GT)
levels(df4$GT) <- c("BY", "het", "RM")

#######################################################################################################
# Proportion of each GT in each clone
ggplot(data = df4, aes(ID, prop, fill = GT)) + 
  geom_bar(stat = "identity", position = "dodge") + facet_grid(algn~.)

#######################################################################################################
# Absolute GT counts for each clone
ggplot(data = df4, aes(ID, count, fill = GT)) + 
  geom_bar(stat = "identity", position = position_dodge2(width = 0.5, preserve = "single")) + 
  scale_y_log10() + labs(color="Genotype") + xlab("Ancestral Clone ID") +
  facet_grid(algn~.) +
  ggtitle("Genotype counts for each ancestral clone when reads are aligned to each reference sequence")

#######################################################################################################
# Read depth frequencies comparing the refseq alignments among the ancestral clones
ggplot(anc_all_hiQual) +
  geom_histogram(aes(x=Sum_DP, color=algn), position="identity", bins=50, alpha=0) + 
  facet_wrap(~ID, ncol=4) + xlim(c(0,250)) + labs(color="Alignment") + xlab("Read Depth")

#######################################################################################################


anc_all_hiQual_het <- anc_all_hiQual %>% filter(GT=="0/1")

#######################################################################################################
# Plot total reads supporting BY and RM alleles for each clone and each alignment
df5 <- anc_all_hiQual_het %>% 
  dplyr::group_by(algn, ID) %>% 
  dplyr::summarize(sum_BY=sum(Ref_DP, na.rm=T), sum_RM=sum(Alt_DP, na.rm=T))

mypalette<-c(brewer.pal(8,"Dark2"), brewer.pal(4,"Set1"), brewer.pal(8,"Set2"))

df5 %>% 
  filter(ID != "N_G00") %>%
  ggplot() +
    geom_point(aes(x=sum_BY, y=sum_RM, color=ID, shape=algn), size = 3, alpha=0.8) + 
    geom_abline(aes(slope=1, intercept=0), alpha=0.7) + 
    labs(shape="Alignment") + xlab("# reads supporting BY allele") + ylab("# reads supporting RM allele") + 
    scale_colour_manual(values=mypalette) +
    xlim(c(0,1E6)) + ylim(c(0,1E6)) #+
    #ggtitle("Total number of reads in each ancestral clone supporting each allele for sites called heterozygous when reads are aligned to each reference sequence")

#######################################################################################################
# Distribution of proportion of reads supporting BY in heterozygous calls
anc_all_hiQual_het %>% 
  # filter(ID != "F_B00" & ID != "F_E00") %>%
ggplot() +
  geom_vline(aes(xintercept=0.5), color="gray50") +
  geom_freqpoly(aes(x=fract_Ref, color=algn, fill=algn), position="identity", bins=21) + 
  theme_minimal() +
  facet_grid(Line~Tx) +
  # facet_wrap(~ID, ncol=6) + 
  scale_color_manual(values=c("brown", "gray20", "blue")) +
  scale_x_continuous(breaks=c(0, 0.25, 0.5, 0.75, 1)) + 
  xlab("Fraction of reads supporting BY")

anc_het_fract_BY <- anc_all_hiQual_het %>% 
  # filter(ID != "F_B00" & ID != "F_E00") %>% 
  group_by(algn, ID) %>%
  summarize(mean_fract_BY = mean(fract_Ref, na.rm=T), stdDev = sd(fract_Ref, na.rm=T))

anc_het_fract_BY %>%
  filter(ID != "F_B00" & ID != "F_E00") %>%
ggplot() + geom_boxplot(aes(algn, mean_fract_BY))

# Compute the analysis of variance
algn_aov <- anc_het_fract_BY %>%
  filter(ID != "F_B00" & ID != "F_E00") %>% aov(mean_fract_BY ~ algn, data=.)

summary(algn_aov)

anc_BYm_fract_BY <- anc_het_fract_BY %>%
  filter(algn == "BYm" & ID != "F_B00" & ID != "F_E00")

t.test(x=anc_BYm_fract_BY$mean_fract_BY, alternative = "two.sided",
       mu = 0.5,conf.level = 0.95)

anc_het_fract_BY %>%
  filter(ID != "F_B00" & ID != "F_E00") %>% group_by(algn) %>% summarize(meanofmeans = mean(mean_fract_BY))


#######################################################################################################
# Find consensus sites

anc_BYm_hiQual <- anc12_BYm_hiQual
anc_BYm_hiQual$GT <- droplevels(anc_BYm_hiQual$GT)
anc_BYm_hiQual$GTx <- as.numeric(anc_BYm_hiQual$GT)-2

anc_het_wide <- anc_BYm_hiQual %>% filter(GT == "0/0" | GT == "0/1" | GT == "1/1") %>%
  pivot_wider(id_cols=c(CHROM, POS, POSi, REF, Alt1), values_from=GTx, 
              names_from=ID, names_prefix="GT.", values_fill=NA)

str(anc_het_wide)

anc_het_wide$SumGT <- anc_het_wide %>% select(.,ends_with("00")) %>% rowSums(na.rm=T)

anc_het_wide$countNA <- anc_het_wide %>% select(.,ends_with("00")) %>% is.na() %>% rowSums()

anc_het_wide %>% filter(SumGT != 0) %>% View 

anc_sum_noCall <- anc_het_wide %>% select(.,ends_with("00")) %>% is.na() %>% colSums()

anc_het_cnsns <- anc_het_wide %>% filter(SumGT==0 & countNA<5) %>% select(., !ends_with("00"))

##################
# Export position table for HaplotypeCaller

##################
# Found two possible gene conversions, one in F_A00 and N_D00
anc_het_wide %>% filter(POSi >= 3676574 &  POSi <= 3682654) %>% View

anc12_BYm_hiQual %>% 
  filter(Tx=="N") %>% 
  filter(POSi >= 3654000 & POSi <= 3800000) %>%
  ggplot() + geom_jitter(aes(x=POSi, y=GT, color=GT), width=0, height=0.1) +
  facet_grid(ID~.)


anc_BYm_hiQual %>% filter(ID=="N_D00", GTx != 0)

anc_BYm_hiQual <- anc_BYm_hiQual %>% separate(ID, sep="_", c("Tx", "Line"), remove=F)


# Get chromosome boundaries for plots ################################

chrom_lngths <- c(230218, 813184, 316620, 1531933, 576874, 270161, 1090940, 562643, 
                  439888, 745751, 666816, 1078177, 924431, 784333, 1091291, 948066)
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
                          CHROM=c(1:16), 
                          stringsAsFactors = F)
  return(chr_bound)
}
chrom_bound <- chrom_coord(chrom_lngths)

g_lngth <- sum(chrom_lngths)

##

###########################################
# Plot genotypes across genome positions
anc_BYm_RMxBY <- anc_BYm_hiQual %>% 
  filter(GT == "0/0" | GT == "0/1" | GT == "1/1")

t_GTx <- anc12_BYm_hiQual %>% 
  dplyr::count(Tx,Line,GT) %>% 
  dplyr::group_by(Line) %>% 
  dplyr::mutate(prop=prop.table(n))

ggplot() +
  geom_point(data=anc12_BYm_hiQual, aes(x=POSi, y=as.numeric(GT), color=as.factor(GT)), size=1) +
  theme_minimal() + scale_color_manual(values=c("orangered3","black","deepskyblue4","red"), guide=F) +
  theme(axis.text.y = element_text(size=8), axis.text.x = element_text(size=6), 
        panel.grid.major.x=element_blank(), panel.grid.minor.x=element_blank()) +
  scale_fill_manual(values=c("orangered3","black","deepskyblue4","red"), guide=) +
  scale_y_discrete(breaks=c(-1,0,1,2), labels=c("BY","het","RM","GT_freq")) + 
  scale_x_continuous(labels=as.character(chrom_bound$CHROM), breaks=(chrom_bound$Start+chrom_bound$End)/2) +
  #ylab(NA) + 
  theme(axis.title.y = element_blank()) +
  xlab("Genome Position") + labs(title="") +
  guides(col=guide_legend(title="Genotype")) +
  annotate(geom="rect", xmin=chrom_bound$Start[c(TRUE, FALSE)], 
           xmax=chrom_bound$End[c(TRUE, FALSE)], ymin=-Inf, ymax=Inf, alpha=0.2) + 
  facet_grid(Line~Tx, scales = "fixed") + 
  theme(panel.background = element_rect(fill = "white", color="gray70"), 
        strip.background = element_rect(color="gray70", fill="white"))

anc_BYm_RMxBY %>% 
  filter(Tx == "F") %>% 
  filter(Line %in% c("A00", "B00", "C00", "E00")) %>% 
ggplot() +
  geom_point(aes(x=POSi, y=fract_Ref), size=0.2, shape=1) +
  theme_minimal() + # scale_color_manual(values=c("orangered3","black","deepskyblue4","red"), guide=F) +
  theme(axis.text.y = element_text(size=8), axis.text.x = element_text(size=6), 
        panel.grid.major.x=element_blank(), panel.grid.minor.x=element_blank()) +
  scale_x_continuous(labels=as.character(chrom_bound$CHROM), breaks=(chrom_bound$Start+chrom_bound$End)/2) +
  theme(axis.title.y = element_blank()) +
  xlab("Genome Position") + labs(title="") +
  # guides(col=guide_legend(title="Genotype")) +
  annotate(geom="rect", xmin=chrom_bound$Start[c(TRUE, FALSE)], 
           xmax=chrom_bound$End[c(TRUE, FALSE)], ymin=-Inf, ymax=1, alpha=0.2) + 
  facet_grid(.~Line, scales = "fixed") + 
  theme(panel.background = element_rect(fill = "white", color="gray70"), strip.background = element_rect(
    color="gray70", fill="white"))




