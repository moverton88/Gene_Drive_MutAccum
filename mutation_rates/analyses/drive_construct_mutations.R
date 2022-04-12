# Mutations may have occurred in the drive construct, disabling
# its components, and confounding any lack of difference observed
# in mutation rates. Therefore, reads from each end-point clone were 
# aligned with that strain's drive cassette and variants were called.

D_const_variants <- HC_multiVcfParse_allVar(D_drive_const_file, all_Alts = T) 
C_const_variants <- HC_multiVcfParse_allVar(C_drive_const_file, all_Alts = T) 
W_const_variants <- HC_multiVcfParse_allVar(W_drive_const_file, all_Alts = T)


drive_const_variants_fltr <- drive_const_variants %>% 
  filter(Sum_DP >= 6) %>% 
  filter(QUAL >= 100)


alt_NAs <- drive_const_variants %>% 
  # select(contains("ALT", ignore.case = F)) %>% 
  apply(., MARGIN = 2, function(x) sum(is.na(x)))

n_rows <- nrow(drive_const_variants)

cols_remove <- names(alt_NAs)[alt_NAs == n_rows]

drive_const_variants %>% filter(!is.na(ALT6), GT != "0/0")
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

anc_poly <- drive_const_variants %>% filter(GT != "0/0", GT != "./.", Rep == "00") %>% distinct(POSi)
drive_const_variants %>% filter(GT != "0/0", GT != "./.", Rep != "00", !POSi %in% anc_poly$POSi)


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





