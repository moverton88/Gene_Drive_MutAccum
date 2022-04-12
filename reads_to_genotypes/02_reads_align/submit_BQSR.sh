#!/bin/bash

##################################################################
# CHOOSE REFERENCE SEQUENCE ######################################
# RM = RM reference
# BY = BY reference
# BYm = BY masked reference
export alignRef=RM
export alignTag="_local"
# _local

if [ ${alignRef} == BY ]; then
    export bamDir=/oasis/tscc/scratch/mioverto/geneDrive/dualRef/BY_aligned/bam
    export refSeq=/home/mioverto/geneDrive/refseq/BY/S288C_R64_refseq.fna
    export vcfIn=/home/mioverto/geneDrive/POS_files/RMvcf/RMxBY_ref_noMit.vcf
    elif [ ${alignRef} == RM ]; then
    export bamDir=/oasis/tscc/scratch/mioverto/geneDrive/dualRef/RM_aligned/bam
    export refSeq=/home/mioverto/geneDrive/refseq/RM/RM_refseq_UCSD_2020_v4.fna
    export vcfIn=/home/mioverto/geneDrive/POS_files/RMxBY_ref_rev.vcf
    else
    echo "reference does not exist"
fi

export DATE=$(date +'%m_%d_%Y')
export script=/home/mioverto/code/align/calcAndApplyBQSR.sh
export logDir=/oasis/tscc/scratch/mioverto/geneDrive/log/BSQR


for bamfile in ${bamDir}/DeDup/N_B*${alignRef}.dm.bam; do
    # export R1FILE=/oasis/tscc/scratch/mioverto/data/MAseq1/reads/trim/half-L100_1_R1P.trimmed.fastq
    export bamDeDup=${bamfile}
    export tmp=$(basename "${bamDeDup}" .dm.bam)
    export index=${tmp:0:5}_${alignRef}
    export bamTable=${bamDir}/BQSR_tables/${index}_bqsr${alignTag}.table
    export bamCrct=${bamDir}/DeDup_crct/${index}${alignTag}.crct.bam
    export bamOut=${bamDir}/BQSR_bam/${index}_bqsr${alignTag}.bam
    # export VCFout=${VCFdir}/${index}.vcf
    echo "Submitting BSQR ${tmp}"
# done
    qsub \
        -V \
        -N BQSR_${index} \
        -o ${logDir}/${index}_${alignRef}_${DATE}.out \
        -e ${logDir}/${index}_${alignRef}_${DATE}.err \
        ${script}
done
