#!/bin/bash
#PBS -A mioverto
#PBS -l nodes=1,walltime=3:00:00


# Switching to java 1.8
PATH=/usr/lib/jvm/java-1.8.0-openjdk-1.8.0.252.b09-2.el7_8.x86_64/bin:$PATH

# R1COMP=$FQDIR/F_A01_1_R1.fastq.gz
# R2COMP=$FQDIR/F_A01_1_R2.fastq.gz

# Record inputs to log
echo fastq inputs:; echo $R1COMP; echo $R2COMP
echo Illumina adapter: $ADAPTER
echo Trimmomatic location: $TRIMMO

# R1COMP=$readsRawDir/F_A00_1_R1.fastq

if [[ ${R1COMP} =~ .*gz.* ]] 
then

   echo "${R1COMP} is compressed, gunzipping."
   gunzip ${R1COMP}
   R1FILE=${R1COMP/.gz/}
 else
    echo "${R1COMP} is already uncompressed."
    R1FILE=${R1COMP}
fi

if [[ ${R2COMP} =~ .*gz.* ]] 
then
   echo "${R2COMP} is compressed, gunzipping."
   gunzip ${R2COMP}
   R2FILE=${R2COMP/.gz/}
 else
    echo "${R2COMP} is already uncompressed."
    R2FILE=${R2COMP}
fi

# Generate the output file paths for paired and unpaired files
fnmR1=${R1FILE}
tmpR1=${trimDir}/$(basename ${fnmR1})

trimR1P=${tmpR1/R1.fastq/R1P.trim.fastq}
trimR2P=${tmpR1/R1.fastq/R2P.trim.fastq}
trimR1U=${tmpR1/R1.fastq/R1U.trim.fastq}
trimR2U=${tmpR1/R1.fastq/R2U.trim.fastq}
# echo $trimR1P

# Send Trimmomatic's logs to same folder; changing file ext later
tmpLog=${logDir}/$(basename ${fnmR1})
trimLog=${tmpLog/R1.fastq/trimlog.txt}



java -jar ${TRIMMO} PE \
   -trimlog ${trimLog} \
   $R1FILE $R2FILE \
   ${trimR1P} ${trimR1U} \
   ${trimR2P} ${trimR2U} \
   ILLUMINACLIP:${ADAPTER}:2:30:8 \
   LEADING:10 \
   TRAILING:10 \
   HEADCROP:15 \
   SLIDINGWINDOW:6:10 \
   MINLEN:30

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