##

repeatBed <- function(tbl_flnm) {
  require(bedr)
  RMxBYvcf <- read.vcf(RMxBY_flnm)
  RMxBYtbl <- vcf2bed(RMxBYvcf)
  rptTbl <- read.table(tbl_flnm, header=F, skip=3)[, 5:7]
  zeroStartPos <- rptTbl[,2] - 1
  rptTbl[,2] <- zeroStartPos
  colnames(rptTbl) <- colnames(RMxBYtbl)
  BedTbl <- rbind(rptTbl, RMxBYtbl)
  BedTbl_srt <- BedTbl[with(BedTbl, order(chr, start)), ]
  return(BedTbl_srt)
  write.table(BedTbl_srt, file = rptBedOut, sep="\t", col.names=F, row.names = F, quote=F)
}

tableToIntervals <- function(rTable) {
  require(bedr)
  # RMxBYvcf <- read.vcf(RMxBY_flnm)
  # RMxBYtbl <- vcf2bed(RMxBYvcf)
  # rptTbl <- read.table(tbl_flnm, header=F, skip=3)[, 5:7]
  zeroStartPos <- rTable$POSi - 1
  rTable$POSi <- zeroStartPos
  colnames(rptTbl) <- colnames(RMxBYtbl)
  BedTbl <- rbind(rptTbl, RMxBYtbl)
  BedTbl_srt <- BedTbl[with(BedTbl, order(chr, start)), ]
  return(BedTbl_srt)
  # write.table(BedTbl_srt, file = rptBedOut, sep="\t", col.names=F, row.names = F, quote=F)
}