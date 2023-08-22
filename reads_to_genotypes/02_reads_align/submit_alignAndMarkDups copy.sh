#!/bin/bash

##################################################################
# CHOOSE SEQUENCING RUN DATASET ######################################
# MAseq1
# MAseq2
# MAseq3
# export seqRun=MAseq1

##################################################################
# CHOOSE REFERENCE SEQUENCE ######################################
# RM = RM reference
# BY = BY reference
# BYm = BY reference with RM variant sites masked with Ns
# Cas9 = Cas9 gRNA construct
# W303 = 
# YJM = YJM
export ref=W303

export proj=Sui_etal_2020
 
export readsDir=/oasis/tscc/scratch/mioverto/LOH_methods/${proj}/reads/trim
# export readsDir=/oasis/tscc/scratch/mioverto/geneDrive/reads/trim

bamDirPre=/oasis/tscc/scratch/mioverto/LOH_methods/${proj}
# bamDirPre=/oasis/tscc/scratch/mioverto/geneDrive/dualRef/
# bamDirPre=/oasis/tscc/scratch/mioverto/geneDrive/Cas9/bam
if [ ${ref} == RM ]; then
    export refSeq=/home/mioverto/geneDrive/refseq/RM/RM_refseq_UCSD_2020_v4.fna
    export bamDir=${bamDirPre}/${ref}_aligned/bam
    elif [ ${ref} == BY ]; then
    export refSeq=/home/mioverto/geneDrive/refseq/BY/S288C_R64_refseq.fna
    export bamDir=${bamDirPre}/${ref}_aligned/bam
    elif [ ${ref} == W303 ]; then
    export refSeq=/home/mioverto/LOH_methods/refseq/W303/W303_refseq.fna
    export bamDir=${bamDirPre}/${ref}_aligned/bam
    elif [ ${ref} == YJM ]; then
    export refSeq=/home/mioverto/LOH_methods/refseq/YJM789/YJM789_refseq.fna
    export bamDir=${bamDirPre}/${ref}_aligned/bam
    elif [ ${ref} == BYm ]; then
    export refSeq=/home/mioverto/geneDrive/refseq/BYm/S288C_R64_masked.fna
    export bamDir=${bamDirPre}
    elif [ ${ref} == Cas9_N ]; then
    export refSeq=/home/mioverto/geneDrive/refseq/Drive/GFP_NrsR_cassette.fna
    export bamDir=${bamDirPre}
    elif [ ${ref} == Cas9_H ]; then
    export refSeq=/home/mioverto/geneDrive/refseq/Drive/Cas9_GFP_NrsR.fna 
    export bamDir=${bamDirPre}
    elif [ ${ref} == Cas9_F ]; then
    export refSeq=/home/mioverto/geneDrive/refseq/Drive/gRNA_Cas9_GFP_NrsR_WTade2.fna
    export bamDir=${bamDirPre}
    else
    echo "reference does not exist"
fi

export refPrefix=${REFSEQ/.fna/}

# export script=/home/mioverto/code/align/alignAndMarkDups_vs.sh
export script=/home/mioverto/code/align/alignAndMarkDups.sh
export logDir=/oasis/tscc/scratch/mioverto/geneDrive/log
export DATE=$(date +'%m_%d_%Y')

export metrics=${bamDir}/picard_metrics

# [ -d "${bamDir}/picard_metrics" ] && echo "Directory /path/to/dir exists."


# export R1FILE=${readsDir}/H_A00_1_R1P.trim.fastq
# Submitting jobs in a loop for files that have not been created yet

for R1FILE in ${readsDir}/L*R1P.trim.fastq; do
    # export R1FILE=/oasis/tscc/scratch/mioverto/data/MAseq1/reads/trim/half-L100_1_R1P.trim.fastq
    export R1PFILE=${R1FILE}
    export R1UFILE=${R1FILE/R1P/R1U}
    export R2PFILE=${R1FILE/R1P/R2P}
    export R2UFILE=${R1FILE/R1P/R2U}
    export tmp=$(basename "${R1FILE}" .trim.fastq)
    export index=${tmp%%_*}_${ref}
    # export index=${tmp:0:5}_${ref}
    export bamRaw=${bamDir}/raw/${index}_local.bam
    export bamDeDup=${bamDir}/DeDup/${index}_local.dm.bam
    echo "Submitting ${index}"
# done
    qsub \
        -V \
        -N align_${index} \
        -o ${logDir}/align-bam_${index}_${DATE}.out \
        -e ${logDir}/align-bam_${index}_${DATE}.err \
        ${script}
done

