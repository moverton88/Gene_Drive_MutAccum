#!/bin/bash

##################################################################
# CHOOSE REFERENCE SEQUENCE ######################################
# RM = RM reference
# BY = BY reference
# BYm = BY masked reference
# Cas9 = gRNA and Cas9 integrated construct
export alignRef=Cas9_F
export callRef=${alignRef}

export tag=default 
# export alignTag=""
# _local

if [ ${alignRef} == BY ] ; then
export refSeq=/home/mioverto/geneDrive/refseq/BY/S288C_R64_refseq.fna
export bamDir=/oasis/tscc/scratch/mioverto/geneDrive/dualRef/${alignRef}_aligned/bam
        export gVCFdir=/oasis/tscc/scratch/mioverto/geneDrive/dualRef/${alignRef}_aligned/variants/gVCFs/${callRef}_call
    elif [[ ${alignRef} == RM ]] ; then
    export refSeq=/home/mioverto/geneDrive/refseq/RM/RM_refseq_UCSD_2020_v4.fna
        export bamDir=/oasis/tscc/scratch/mioverto/geneDrive/dualRef/${alignRef}_aligned/bam
        export gVCFdir=/oasis/tscc/scratch/mioverto/geneDrive/dualRef/${alignRef}_aligned/variants/gVCFs/${callRef}_call
    elif [[ ${alignRef} == Cas9_N ]] ; then
        export bamDir=/oasis/tscc/scratch/mioverto/geneDrive/Cas9/bam
        export gVCFdir=/oasis/tscc/scratch/mioverto/geneDrive/Cas9/variants/gVCFs/
        export refSeq=/home/mioverto/geneDrive/refseq/Drive/GFP_NrsR_cassette.fna
    elif [[ ${alignRef} == Cas9_H ]] ; then
        export bamDir=/oasis/tscc/scratch/mioverto/geneDrive/Cas9/bam
        export gVCFdir=/oasis/tscc/scratch/mioverto/geneDrive/Cas9/variants/gVCFs/
        export refSeq=/home/mioverto/geneDrive/refseq/Drive/Cas9_GFP_NrsR.fna
    elif [[ ${alignRef} == Cas9_F ]] ; then
        export bamDir=/oasis/tscc/scratch/mioverto/geneDrive/Cas9/bam/
        export gVCFdir=/oasis/tscc/scratch/mioverto/geneDrive/Cas9/variants/gVCFs/
        export refSeq=/home/mioverto/geneDrive/refseq/Drive/gRNA_Cas9_GFP_NrsR_WTade2.fna
    else
        echo "alignment reference does not exist"
fi

```
module load samtools
PATH=/usr/lib/jvm/java-1.8.0-openjdk-1.8.0.252.b09-2.el7_8.x86_64/bin:$PATH
GATK=/home/mioverto/bin/gatk_4/gatk-package-4.1.8.0-local.jar

samtools faidx $refSeq
refdict=${refSeq/.fna/.dict}
# GATK requires its own indexing prep
java -jar $GATK CreateSequenceDictionary \
    -R $refSeq \
    -O $refdict

```
export metricsDir=${bamDir}/depth_metrics/

export DATE=$(date +'%m_%d_%Y')
export script=/home/mioverto/code/variants/01_call_gVCF_gatkHC_default.sh
export logDir=/oasis/tscc/scratch/mioverto/geneDrive/log/callGvcf
# export vcfIn=/home/mioverto/geneDrive/POS_files/RMvcf/RMxBY_ref_noMit.vcf
# export vcfIn=/home/mioverto/geneDrive/POS_files/RMxBY_ref_rev.vcf

for bamfile in ${bamDir}/DeDup/F*_local.dm.bam; do
    # export R1FILE=/oasis/tscc/scratch/mioverto/data/MAseq1/reads/trim/half-L100_1_R1P.trimmed.fastq
    export bamDeDup=${bamfile}
    export tmp=$(basename "${bamDeDup}" .dm.bam)
    export index=${tmp:0:5}_${callRef}
    export bamAlgnMetrics=${metricsDir}/${index}_${tag}.algn.txt
    export bamWGSmetrics=${metricsDir}/${index}_${tag}.wgs.txt
    export gVCFout=${gVCFdir}/${index}_${tag}.g.vcf
    # export VCFout=${VCFdir}/${index}.vcf
    echo "Submitting call gVCF ${tmp}"
# done
    qsub \
        -V \
        -N call_${index} \
        -o ${logDir}/gVCF_${index}_${DATE}.out \
        -e ${logDir}/gVCF_${index}_${DATE}.err \
        ${script}
done

