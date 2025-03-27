#!/bin/bash
#PBS -A mioverto
#PBS -l nodes=1 
#PBS -l walltime=5:00:00

###############################################################################
# 02_trim_reads.sh
# Takes in a set of either single- or paired-end reads and removes erroneous 
# bases from each end of the read. lead and trail remove bases from each end
# based on quality, head removes the given number of bases from the first 
# position in the read, crop removes bases from the end until the read is 
# of the given legnth. An optional adapter file can be used to remove known 
# adapter sequences from the reads as well.
#*****************************************************************************#


while getopts 'r:m:o:C:L:T:H:A:g:' OPTION; do
    case "$OPTION" in
    r)
        reads1="$OPTARG"
        echo "Pair 1 Reads file $OPTARG"
    ;;

    m)
        mode="$OPTARG"
        if [[ -z ${mode} ]]; then
            printf "Defaulting to Paired-end mode"
            mode=Paired
        else
            echo "${mode} mode selected"
        fi
    ;;  

    o)
        outputDir="$OPTARG"
        echo "The output directory is $OPTARG"
    ;;   

    C)
        crop="$OPTARG"
         if [[ -z ${crop} ]]; then
            printf "Defaulting to 1000bp max read length"
            crop=1000
        else
            echo "Max read length ${crop}"
        fi
    ;; 

    L)
        lead="$OPTARG"
        echo "Minimum quality for inclusion at read start: $OPTARG"
    ;;    

    T)
        trail="$OPTARG"
        echo "Minimum quality for inclusion at read end: $OPTARG "
    ;;    

    H)
        head="$OPTARG"
        echo "Crop $OPTARG bases from the start of the read"
    ;;    


    A)
        adapters="$OPTARG"
        echo "Adapter file $OPTARG"
    ;;

    ?)
        echo "script usage: $(basename \$0) [-r reads1File] [-m modeSingleOrPaired] [-o outputDirectory]
        [-L leadingMinQual] [-T TrailingMinQual] [-H headCrop] [-C maxReadLength] [-A adapterSeqFile] 
        [-g logDirectory]"
    esac
done

if [[ ${mode} == "Paired" || ${mode} == "paired" || ${mode} == "p" ]]; then
    reads2=${reads1/R1/R2}
    trimOpt=PE
else
    trimOpt=SE
fi

#  adapters=/home/mioverto/bin/Trimmomatic-0.36/adapters/NexteraPE-PE.fa

# Record inputs to log
echo fastq inputs:; echo $reads1; echo $reads2
echo Illumina adapter: $adapters
echo Trimmomatic location: $trimApp

if [[ ${reads1} =~ .*gz.* ]] 
then
   echo "${reads1} is compressed, gunzipping."
   gunzip ${reads1}
   R1file=${reads1/.gz/}
 else
    echo "${reads1} is already uncompressed."
    R1file=${reads1}
fi

if [[ ${reads2} =~ .*gz.* ]] 
then
   echo "${reads2} is compressed, gunzipping."
   gunzip ${reads2}
   R2file=${reads2/.gz/}
 else
    echo "${reads2} is already uncompressed."
    R2file=${reads2}
fi

# Generate the output file paths for paired and unpaired files
tmpReadFile=${outputDir}/$(basename ${R1file})

trimR1=${tmpReadFile/R1.fastq/R1.trim.fastq}
trimR2=${tmpReadFile/R1.fastq/R2.trim.fastq}
trimR1U=${tmpReadFile/R1.fastq/R1U.trim.fastq}
trimR2U=${tmpReadFile/R1.fastq/R2U.trim.fastq}
# echo $trimR1

# Send Trimmomatic's logs to same folder; changing file ext later
trimLog=${readsDir}/trim_log/$(basename ${R1file/_R1.fastq/})_trimlog.txt

if [ ! -f ${adapters} ]; then
    java -jar ${trimApp} ${trimOpt} \
        -trimlog ${trimLog} \
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
        -trimlog ${trimLog} \
        $R1file $R2file \
        ${trimR1} ${trimR1U} \
        ${trimR2} ${trimR2U} \
        ILLUMINACLIP:${adapters}:2:30:8 \
        LEADING:${lead} \
        TRAILING:${trail} \
        HEADCROP:${head} \
        CROP:${crop} \
        SLIDINGWINDOW:6:10 \
        MINLEN:30
fi


#  ILLUMINACLIP:${ADAPTER}:1:30:10 \
#     LEADING:3 \
#     TRAILING:3 \
#     SLIDINGWINDOW:10:7 \
#     HEADCROP:10 \
#     CROP:80 \
#     MINLEN:20

# ILLUMINACLIP:TruSeq3-PE.fa:2:30:10
#     LEADING:3 
#     TRAILING:3 
#     SLIDINGWINDOW:4:15 
#     MINLEN:36