#!/bin/bash
 
# This script has seven requirements:
# 1 - A directory where all read files are stored
# 2 - The expected length of the reads
# 3 - Paired read files that have the same sample_ID and are distinguished by 
#   R1 and R2
# 4- A reference sequence file with leading "Name_" 
#   and a .fasta suffix in the filename (eg. path/BY_refseq.fasta)
# 5 - A directory for the alingment .bam files
# 6 - The alignment script path and filename. 
#   I keep all my scripts for each project in a directory on my home drive
# 7 - A path to a log directory to store stdout and stderr outputs

# Submitting jobs in a loop for files that have not been created yet

#!!! Change to your account username
export usrname=mioverto

# replace path/to/... with the actual path to each
export readsDir=/oasis/tscc/scratch/${usrname}/path/to/reads
# export readsDir=/oasis/tscc/scratch/mioverto/LOH_methods/Overton_etal/reads

export alignDir=/oasis/tscc/scratch/${usrname}/path/to/alignments
# export alignDir=/oasis/tscc/scratch/mioverto/LOH_methods/Overton_etal/alignments

export refSeq=/oasis/tscc/scratch/${usrname}/path/to/refseq/RefName_refseq.fasta
# export refSeq=/oasis/tscc/scratch/mioverto/LOH_methods/Overton_etal/references/final/S288C/S288C_R64_refseq.fasta

export alignScript=/home/${usrname}/project/code/align_markDups.sh
# chmod +x $alignScript
# export alignScript=/home/mioverto/code/align_markDups.sh

export logDir=/oasis/tscc/scratch/${usrname}/path/to/logs
# export logDir=/oasis/tscc/scratch/mioverto/LOH_methods/Overton_etal/logs

export readLength=100

# you can also narrow the operation to a particular set_ID
# or leave empty to operate on all pairs of read files
export set_ID="Test"

for reads1 in ${readsDir}/${set_ID}*R1*fastq; do
    export reads1file=${reads1}
    # export reads1file=${readsDir}/N_C12_R1.trim.fastq
    # Get a sample ID from the file name
    export index=$(basename "${reads1file}" .fastq)
    echo "Submitting align_${index}"
# done
    qsub -V -N align_${index} \
                    -A ${usrname} \
                    -l walltime=08:00:00 \
                    -o ${logDir}/align_${index}_$(date +'%Y_%m_%d').out \
                    -e ${logDir}/align_${index}_$(date +'%Y_%m_%d').err \
                    -F "-r ${reads1file} -a ${alignDir} -R ${refSeq} -L ${readLength}" \
                    ${alignScript}
done

# or run this in interactive mode

${alignScript} -r ${reads1file} -a ${alignDir} -R ${refSeq} -L ${readLength}