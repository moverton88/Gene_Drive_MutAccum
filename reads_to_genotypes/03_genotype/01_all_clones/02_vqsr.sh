#!/bin/bash
#PBS -A mioverto
#PBS -l nodes=1 
#PBS -l walltime=5:00:00

# Combine VCF files of individual clones into one master VCF

PATH=/usr/lib/jvm/java-1.8.0-openjdk-1.8.0.252.b09-2.el7_8.x86_64/bin:$PATH
GATK=/home/mioverto/bin/gatk_4/gatk-package-4.1.8.0-local.jar

alignRef=BY
export refSeq=/home/mioverto/geneDrive/refseq/BY/S288C_R64_refseq.fna
export VCFdir=/oasis/tscc/scratch/mioverto/geneDrive/dualRef/${alignRef}_aligned/variants/allVarVcfs/
export VCFin=${VCFdir}/A_BYaBYc_allVar.vcf
export VCFvqsr=${VCFdir}/A_BYaBYc_allVar_vqsr.vcf
export VCFout=${VCFdir}/A_BYaBYc_allVar_filter.vcf
export VCFrmby=/home/mioverto/geneDrive/POS_files/RMxBY_ref_bcf.vcf
# java -jar $GATK IndexFeatureFile -I ${VCFrmby}
export VCFanc=${VCFdir}/anc_pool_${alignRef}c_filtered.vcf
export VQSRdir=${VCFdir}/vqsr/

java -jar $GATK VariantRecalibrator \
    -R ${refSeq} \
    -V ${VCFin} \
    --resource:RMxBY,known=false,training=true,truth=false,prior=12.0 ${VCFrmby} \
    --resource:Anc,known=false,training=true,truth=true,prior=15.0 ${VCFanc} \
    -an QD -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
    -mode SNP \
    -O ${VQSRdir}/output.recal \
    --tranches-file ${VQSRdir}/output.tranches \
    --rscript-file ${VQSRdir}/output.plots.R

java -jar $GATK ApplyVQSR \
   -R ${refSeq} \
   -V ${VCFin} \
   -O ${VCFvqsr} \
   --truth-sensitivity-filter-level 99.0 \
   --tranches-file ${VQSRdir}/output.tranches \
   --recal-file ${VQSRdir}/output.recal \
   -mode SNP

   java -jar $GATK SelectVariants \
     -R ${refSeq} \
     -V ${VCFin} \
     --exclude-filtered \
     --exclude-non-variants \
     --remove-unused-alternates \
     -O ${VCFout}
