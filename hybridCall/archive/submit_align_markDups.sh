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

export alignDir=/oasis/tscc/scratch/${usrname}/path/to/alignments

export refSeq=/oasis/tscc/scratch/${usrname}/path/to/refseq/RefName_refseq.fasta

export alignScript=/home/${usrname}/project/code/align_markDups.sh
# before script will run, must have permissions set
# run the following once after uploading script file
# chmod +x $alignScript

export logDir=/oasis/tscc/scratch/${usrname}/path/to/logs

export readLength=100

# you can also narrow the operation to a particular set_ID
# or leave empty to operate on all pairs of read files
export set_ID="test"

for reads1 in ${readsDir}/${set_ID}*R1*fastq; do
    export reads1file=${reads1}
    # Get a sample_ID from the file name
    export sample_ID=$(basename "${reads1file}" .fastq)
    echo "Submitting align_${sample_ID}"
# done
    qsub -V -N align_${sample_ID} \
                    -A ${usrname} \
                    -l walltime=08:00:00 \
                    -o ${logDir}/align_${sample_ID}_$(date +'%Y_%m_%d').out \
                    -e ${logDir}/align_${sample_ID}_$(date +'%Y_%m_%d').err \
                    -F "-r ${reads1file} -a ${alignDir} -R ${refSeq} -L ${readLength}" \
                    ${alignScript}
done

# or run this in interactive mode

${alignScript} -r ${reads1file} -a ${alignDir} -R ${refSeq} -L ${readLength}