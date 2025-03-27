#!/bin/bash
#SBATCH --partition=hotel
#SBATCH --qos=hotel
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=12:00:00
#SBATCH --account=htl100
#SBATCH --get-user-env
#SBATCH --job-name=trim_reads


###############################################################################
# 02_trim_reads.sh
# Takes in a set of either single- or paired-end reads and removes erroneous 
# bases from each end of the read. lead and trail remove bases from each end
# based on quality, head removes the given number of bases from the first 
# position in the read, crop removes bases from the end until the read is 
# of the given legnth. An optional adapter file can be used to remove known 
# adapter sequences from the reads as well.
#*****************************************************************************#


# reads1=${readsDir}/raw/F_A00_R1.fastq.gz  
# trimOpt=PE
# outputDir=${readsDir}
# crop=1000
# lead=10
# trail=10
# head=10
# adapters=/tscc/nfs/home/mioverto/bin/Trimmomatic-0.39/adapters/NexteraPE-PE.fa



# Generate the output file paths for paired and unpaired files
R1file=${reads1}
R2file=${reads1/_R1/_R2}
readsFile=$(basename ${reads1})
tmpReadFile=${outputDir}/${readsFile}

trimR1=${tmpReadFile/R1.fastq/R1.fastq}
trimR2=${tmpReadFile/R1.fastq/R2.fastq}
trimR1U=${tmpReadFile/R1.fastq/R1U.fastq}
trimR2U=${tmpReadFile/R1.fastq/R2U.fastq}
# echo $trimR1

# Send Trimmomatic's logs to same folder; changing file ext later
# trimLog=${readsDir}/trim_log/${readsFile%%.*}_trimlog.txt

if [ ! -f ${adapters} ]; then
    java -jar ${trimApp} ${trimOpt} \
        $R1file $R2file \
        ${trimR1} ${trimR1U} \
        ${trimR2} ${trimR2U} \
        LEADING:${lead} \
        TRAILING:${trail} \
        HEADCROP:${head} \
        CROP:${crop} \
        SLIDINGWINDOW:6:10 \
        MINLEN:30
else
    java -jar ${trimApp} ${trimOpt} \
        $R1file $R2file \
        ${trimR1} ${trimR1U} \
        ${trimR2} ${trimR2U} \
        ILLUMINACLIP:${adapters}:2:30:10 \
        LEADING:${lead} \
        TRAILING:${trail} \
        HEADCROP:${head} \
        CROP:${crop} \
        SLIDINGWINDOW:6:10 \
        MINLEN:30
fi
