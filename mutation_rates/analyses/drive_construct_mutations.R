# Mutations may have occurred in the drive construct, disabling
# its components, and confounding any lack of difference observed
# in mutation rates. Therefore, reads from each end-point clone were 
# aligned with that strain's drive cassette and variants were called.

feature_names <- c("pSNR52", "Target sequence", "gRNA", "tSUP4", 
                   "pTEF1", "Cas9", "tCYC1", "pRPL39", "GFP", "tCYC1",
                   "pTEF1", "NsrR", "tTEF")
feature_labels <- c("", "", "gRNA", "", 
                   "", "Cas9", "", "", "GFP", "",
                   "", "NsrR", "")
feature_type <- c("p", "g", "g", "t", "p", "g", "t", "p", "g", "t", "p", "g", "t")
feature_starts <- c(1, 270, 290, 369, 399, 829, 4981, 5230, 5830, 6549, 6821, 7224, 7802)
feature_ends <- c(269, 289, 368, 388, 797, 4932, 5228, 5823, 6546, 6798, 7219, 7796, 7999)
const_features <- data.frame(feature = feature_names, f_type = feature_type, 
                             start = feature_starts, end = feature_ends)
const_features$feature <- factor(const_features$feature)
const_features$f_type <- factor(const_features$f_type, 
                                levels = c("p", "g", "t"), 
                                labels = c("Constitutive Promoter", "Coding Region", "Terminator"))

const_features$f_label <- feature_labels
const_features$mid <- (const_features$start + const_features$end)/2
const_features$label.y <- ifelse(const_features$end - const_features$start < 350, 1.5, 1.5)
const_features$label.yadj <- c(0.5, -0.5, 0, -0.5, 0.5, 1.5, 0.25, -0.75, 1.5, 0.25, -0.75, 1.5, -0.25)
# const_features$label.yadj <- const_features$label.y + c(-1, -2, -3, -2, rep(0, nrow(const_features) - 3))
const_features$label.xadj <- const_features$mid + c(-450, -350, 0, 400, 600, rep(0, nrow(const_features) - 5))

const_features$segment.y <- const_features$label.yadj + c(0, 0, 0, 0, 0, 0, 0.5, 0.5, 0, 0.5, 0.5, 0, 0.5)
const_features$segment.yend <- 1.5
const_features$segment.x <- const_features$mid + c(0, 0, 0, -1500, -2000, 0, -500, -500, 0, -500, -500, 0, -500)
const_features$segment.xend <- const_features$label.xadj
const_features$fill_color <- c("grey80", "grey60", "grey40")[
                                        as.numeric(const_features$f_type)]
const_features$color_type <- c("Constitutive Promoter", "Coding Region", "Terminator")[
  as.numeric(const_features$f_type)]

D_const_variants <- HC_multiVcfParse_allVar(D_drive_const_file, all_Alts = T)
D_const_variants <- D_const_variants %>% select(!c(ALT3:ALT4, Alt3_DP:Alt6_DP))

C_const_variants <- HC_multiVcfParse_allVar(C_drive_const_file, all_Alts = T)
C_const_variants <- C_const_variants %>% select(!c(Alt3_DP:Alt6_DP))
C_const_variants$POSi <- C_const_variants$POSi + 398
W_const_variants <- HC_multiVcfParse_allVar(W_drive_const_file, all_Alts = T)
W_const_variants <- W_const_variants %>% select(!c(ALT3:ALT4, Alt3_DP:Alt6_DP))
W_const_variants$POSi <- W_const_variants$POSi + 5229

drive_const_variants <- rbind(D_const_variants, C_const_variants, W_const_variants)
drive_const_variants$Tx_ID <- Recode_Tx_ID(drive_const_variants$Tx)
drive_const_variants <- drive_const_variants %>% 
  filter(Sum_DP >= 6, QUAL >= 100, GQ >= 30)


alt_NAs <- drive_const_variants %>% 
  # select(contains("ALT", ignore.case = F)) %>% 
  apply(., MARGIN = 2, function(x) sum(is.na(x)))

n_rows <- nrow(drive_const_variants)

cols_remove <- names(alt_NAs)[alt_NAs == n_rows]

drive_const_variants %>% filter(!is.na(ALT2), GT != "0/0") %>% count(POSi)

drive_const_variants %>% filter(Tx == "F", POSi == 270, !GT %in% c("0/0", "./."))

drive_const_variants %>% 
  filter(POSi > 269, POSi < 291, Tx == "F", GT != "./.") %>% count(POSi, GT)

var_POSi <- drive_const_variants %>% distinct(POSi) %>% pull(POSi)
drive_const_variants %>% 
  # filter(POSi == var_POSi[22]) %>% 
  filter(POSi %in% var_POSi[27:32]) %>% 
  # filter(GT == "1/1")
  filter(!GT %in% c("0/0", "./."))
  # filter(Rep == "00") #%>%
  count(POSi, GT)

anc_poly <- drive_const_variants %>% 
  filter(GT != "0/0", GT != "./.", Rep == "00") %>% 
  distinct(POSi)

const_mut_POSi <- drive_const_variants %>% 
  filter(GT != "0/0", GT != "./.", Rep != "00", !POSi %in% anc_poly$POSi) %>%
  count(POSi) %>% filter(n == 1) %>% pull(POSi)


const_mutants <- drive_const_variants %>% 
  filter(GT != "0/0", GT != "./.", Rep != "00", POSi %in% const_mut_POSi)

const_mutants$type <- factor(ifelse(nchar(const_mutants$REF) > 1 |
                               nchar(const_mutants$ALT1) > 1,
                             "Indel", "SNP"), levels = c("SNP", "Indel"))

const_mutants <- const_mutants %>% 
  mutate(.y = 2.6 + (3 - as.numeric(Tx_ID)),
         l = paste0(REF, " > ", ALT1, " ", c("", "", "Non-Syn ", "", "")))


Tx_const <- data.frame(Tx_ID = Recode_Tx_ID(Tx_ID_levels, "Tx_ID"), 
                       start = c(5229, 398, 1), end = 7999, 
                       y.min = c(4.2, 3.2, 2.2), y.max = c(5.0, 4.0, 3.0))


const_mutant_plot <- const_mutants %>% 
  ggplot() + 
  geom_segment(aes(x = -50, xend = 8050, y = 1.5, yend = 1.5), color = "grey10", size = 2) + 
  geom_rect(data = const_features, 
            aes(xmin = start, xmax = end, ymin = 1, ymax = 2, fill = f_type), 
            color = "grey30", size = 0.25) +
  geom_rect(data = Tx_const, aes(xmin = start, xmax = end, ymin = y.min, ymax = y.max), fill = txPal) +
  geom_text(data = Tx_const, aes(x = end + 150, y = (y.min + y.max)/2, label = Tx_ID), size = 9) +
  # geom_segment(data = const_features, aes(x = label.xadj, xend = mid, y = segment.y, yend = segment.yend)) +
  geom_segment(aes(x = 329.0, xend = 329.0, y = 0, yend = 1.5)) +
  geom_label(data = const_features, aes(x = label.xadj, y = label.yadj, label = f_label),
             lineheight = 0.7,
             label.size = 0, size = 8) +
  geom_point(aes(x = POSi, y = .y, shape = type), 
             color = "white", size = 7) +
  geom_label_repel(aes(x = POSi, y = .y, label = l),
                   ylim = c(5.5, 5.5),
                   label.size = 0,
                   size = 7,
                   segment.color = 'grey50') +
  geom_segment(aes(x = 7999 - 1000, xend = 7999, y = 0.4, yend = 0.4), size = 2.5) +
  geom_text(aes(x = 7999 - 175, y = -0), label = "1000bp", size = 5.5) +
  scale_shape_manual(values = c(18, 20), name = "Mutation") +
  scale_fill_manual(values = c("grey85", "grey70", "grey40"), name = "") +
  scale_color_manual(values = txPal, name = "Strain", drop = F, ) +
  guides(shape = "none",
         fill = guide_legend(
           label.position = "left",
           label.hjust = 1)) +
  # ylim(NA, 20) +
  scale_y_continuous(limits = c(NA, 20), expand = c(0.2, 0.2)) +
  scale_x_continuous(expand = c(0.2, 0.2)) +
  theme(panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        legend.text = element_text(size = 24),
        # legend.spacing.x = unit(1.0, 'mm'),
        # legend.key = element_rect(fill = "grey60", color = NA),
        # legend.direction = "horizontal",
        # legend.margin = margin(t = 0, r = 1, b = 0, l = 1, unit = "mm"),
        legend.position = c(0.77, 0.12),
        plot.margin = margin(0, 10, 10, 10))

const_mutant_plot


ggsave(file.path(outIntDir, "construct_mutant_plot_2022_05.png"), 
       plot = const_mutant_plot,
       device = "png",
       width = 20, height = 14,
       units = "in",
       dpi = 600)


# False positive positions:
# 37 - heterozygous in founders
# 122 - heterozygous in founders
# 267 - low quality het calls
# 270 - end of pSNR52 promoter, reads from genomic promoter align here and 
#    yield calls of insertion between promoter and gRNA
# 273, 274 - low quality het calls
# 276, 279 - hom alt in all clones
# 280 - low quality hom alt call
# 282, 285, 288 - hom alt in all clones
# 297 - possible het in F_C10
# 377, 384 - possible het in F_G12
# 399 - Low quality, possible hom alt in Cas9 clones
# 504 - Het in founders, homopolymer region
# 639 - messy calls, homopolymer region
# 650 - two hets, homopolymer region
# 651 - messy calls, homopolymer region
# 802-806 - messy calls, promoter region

# Real Cas9 mutations
# 1462 4bp deletion, H_C02
# 1724 Aln->Thr, SNP, het, H_C12
# 1774 syn, 1782 Ser->Asn, SNPs in Cas9, het, F_E09
# 3164 9bp R-M-K deletion, het, H_A07
# 4617 2bp insertion, het, H_B09





