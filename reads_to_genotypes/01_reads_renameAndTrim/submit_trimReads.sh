#!/bin/bash

# Script for submitting a set of fastq read files and a trimming script to the remote Cluster.
# Set dir variables
# seqRun=MAseq3
export rawDir=/oasis/tscc/scratch/mioverto/geneDrive/reads/raw
#${seqRun}
# readsTrimDir=${readsRawDir/raw/test}
export trimDir=/oasis/tscc/scratch/mioverto/geneDrive/reads/trim
export logDir=/oasis/tscc/scratch/mioverto/geneDrive/log/trim

# Location of the trimmomatic execution dir
export TRIMMO=/home/mioverto/bin/Trimmomatic-0.36/trimmomatic-0.36.jar
# Location of the script that runs trimmomatic
export script=/home/mioverto/code/trim/01_trim_Reads.sh
# Fasta of the Illumina adapter sequence to remove
export ADAPTER=/home/mioverto/bin/Trimmomatic-0.36/adapters/NexteraPE-PE.fa

export DATE=$(date +'%m_%d_%Y')
# r1file=${rawDir}/F_A00_1_R1.fastq.gz

# Index to uniquely name log files.
#i=0
for r1file in ${rawDir}/*R1*; do
    # echo $r1file
# done
    #i=$(($i+1))
    export R1COMP=${r1file}
    export R1COMP=R1P_hello
    export R2COMP=${R1COMP/R1/R2}
    export tmp=$(basename "${R1COMP/_R1/}")
    export index=${tmp:0:5}
    echo submitting ${index} 
# done
    qsub \
        -V \
        -o ${logDir}/trim_${index}_${DATE}.out \
        -e ${logDir}/trim_${index}_${DATE}.err \
        -N trim_${index} \
        ${script}
done



