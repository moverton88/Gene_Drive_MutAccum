#!/bin/bash
#PBS -A mioverto
#PBS -l nodes=1 
#PBS -l walltime=5:00:00

# Combine VCF files of individual clones into one master VCF

PATH=/usr/lib/jvm/java-1.8.0-openjdk-1.8.0.252.b09-2.el7_8.x86_64/bin:$PATH
GATK=/home/mioverto/bin/gatk_4/gatk-package-4.1.8.0-local.jar

export alignRef=RM
# export refSeq=/home/mioverto/geneDrive/refseq/RM/RM_refseq_UCSD_2020_v4.fna
# export refSeq=/home/mioverto/geneDrive/refseq/BY/S288C_R64_refseq.fna

export finalVCFdir=/oasis/tscc/scratch/mioverto/geneDrive/dualRef/${alignRef}_aligned/variants/allVarVcfs
export filterVCFdir=${finalVCFdir}/filtered

export SNPfn=A_${alignRef}a${alignRef}c_SNPs.vcf
export SNPin=${finalVCFdir}/${SNPfn}
export SNPmark=${filterVCFdir}/${SNPfn/.vcf/_marked.vcf}
export SNPfilter=${filterVCFdir}/${SNPfn/.vcf/_filtered.vcf}

export indelfn=A_${alignRef}a${alignRef}c_indel.vcf
export indelIn=${finalVCFdir}/${indelfn}
export indelMark=${filterVCFdir}/${SNPfn/.vcf/_marked.vcf}
export indelFilter=${filterVCFdir}/${indelfn/.vcf/_filtered.vcf}

# SNP filter
   java -jar $GATK VariantFiltration \
      -R $refSeq \
      -V ${SNPin} \
      -O ${SNPmark} \
      --filter-name "QD2" --filter-expression "QD < 2.0" \
      --filter-name "FS60" --filter-expression "FS > 60.0" \
      --filter-name "MQ20" --filter-expression "MQ < 20.0" \
      --filter-name "QUAL30" --filter-expression "QUAL < 30.0" \
      --filter-name "SOR10" --filter-expression "SOR > 10.0" \
      --filter-name "DP6" --filter-expression "DP < 6"

   java -jar $GATK SelectVariants \
      -V ${SNPmark} \
      -O ${SNPfilter} \
      --exclude-filtered true


# Indel filter
   java -jar $GATK VariantFiltration \
      -R $refSeq \
      -V ${indelIn} \
      -O ${indelMark} \
      --filter-name "QD2" --filter-expression "QD < 2.0" \
      --filter-name "FS200" --filter-expression "FS > 200.0" \
      --filter-name "QUAL30" --filter-expression "QUAL < 30.0" \
      --filter-name "ReadPosRankSum" --filter-expression "ReadPosRankSum < -20.0" \
      --filter-name "DP6" --filter-expression "DP < 6"

   java -jar $GATK SelectVariants \
      -V ${indelMark} \
      -O ${indelFilter} \
      --exclude-filtered true

```

--filter-name "QUAL1000" --filter-expression "QUAL < 1000.0" \
--filter-name "MQRS-12.5" --filter-expression "MQRankSum < -12.5" \


Filter

                    SNP     INDEL

FS                  > 60    > 200
ReadPosRankSum      < -8.0  < -20.0
QUAL                < 30.0  < 30.0
SOR                 > 3.0   NONE
MQ                  < 40.0  NONE
MQRankSum           < -12.5 NONE

```