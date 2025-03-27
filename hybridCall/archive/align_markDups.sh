#!/bin/bash

while getopts 'r:a:R:L:' OPTION; do
   case "$OPTION" in
   r)
      reads1="$OPTARG"
      echo "Reads file 1 $OPTARG"
      ;;

   a)
      alignDirOut="$OPTARG"
      echo "alignment directory is $OPTARG"
      ;;

   R)
      refSeq="$OPTARG"
      echo "Reference sequence $OPTARG"
      ;;

    L) 
      readLength=$OPTARG
      echo "Reads length $OPTARG"
    ;;

   ?)
        echo "script usage: $(basename \$0) [-r readsfile1] [-a alignmentDirOut] [-R referenceSequence] [-L readLength]"
    esac
done


module load bwa
module load samtools
module load bcftools
module load bedtools

gatk=/opt/biotools/GenomeAnalysisTK/4.0.11.0/gatk-package-4.0.11.0-local.jar

refPrefix=${refSeq/.fasta/}
refFile=$(basename "${refSeq}")
refName=${refFile%%_*}

# reads1=${reads1file}
reads2=${reads1/_R1/_R2}
# reads2=${reads1/_R1P/_R2P}
readsName=$(basename "${reads1%%_R1*}")

if [ -z "$sampleName" ]; then  
   sampleName=${readsName}
fi

# alignDirOut=${alignDir}/${P1}
alignRaw=${alignDirOut}/${sampleName}.raw.bam
alignDeDup=${alignRaw/.raw.bam/.dedup.bam}
deDupMetrics=${alignDirOut}/metrics/${sampleName}_dupsMetric.txt
depthFile=${deDupMetrics/_dupsMetric.txt/_depth.txt}
tmpDir=${alignDirOut}/${sampleName}_tmp


if [ ! -d ${alignDirOut}/ ]; then
    echo "Create alignment output dir"
    mkdir ${alignDirOut}/
fi

if [ ! -d ${tmpDir} ]; then
    echo "Create alignment temp dir"
    mkdir ${alignDirOut}/${sampleName}_tmp
fi

if [ ! -d ${alignDirOut}/metrics/ ]; then
    echo "Create metrics dir"
    mkdir ${alignDirOut}/metrics/
fi


if [[ ! -f "${refPrefix}.amb" || ! -f "${refPrefix}.bwt" || ! -f "${refPrefix}.sa" ]]; then
   echo "bwa reference index files does not exist. Generating now..."   
   bwa index -p $refPrefix $refSeq

else
   echo "bwa reference index files already exists"

fi

if [[ ! -f "${refSeq}.img" ]]; then
   echo "bwa index image file does not exist. Generating now..."   

   java -jar $gatk BwaMemIndexImageCreator \
     -I ${refPrefix}.fasta \
     -O ${refSeq}.img

else
   echo "bwa index image file already exists"

fi

if [[ ! -f "${refPrefix}.fasta.fai" ]]; then

   samtools faidx $refSeq

else
   echo "Reference index file already exists"
fi

### Align paired and unpaired reads reads with bowtie2. 
# ${sampleName} used for naming alignments. Output is a sorted bam alignment.

# bwa mem -v [verbose level] 
   # align reads with default bwa settings

# samtools view -b
   # convert alignment from .sam to .bam

# samtools sort -o [output file] -T [temp files] -o [outputFile]
   # sort the .bam file according to position

# samtools index [inputFile]
   # Create bam index file

bwa mem -v 2 $refPrefix $reads1 $reads2 | samtools view -b - | \
samtools sort -O bam -T ${tmpDir} -o ${alignRaw} - 
samtools index ${alignRaw}

rm -dr ${alignDirOut}/${sampleName}_tmp

# Remove duplicate reads and create metrics file
java -jar $gatk MarkDuplicates \
   -I ${alignRaw} \
   -O ${alignDeDup} \
   -M ${deDupMetrics} \
   --REMOVE_DUPLICATES true \
   --TMP_DIR ${alignDir}/${refName}/${sampleName}_tmp

# remove raw alignment files
# rm ${alignRaw/.bam/.rg.bam}
# rm ${alignRaw/.bam/.rg.bam.bai}

# index bam alignment
samtools index ${alignDeDup}

# Generate depth metrics file
java -jar -Xmx16g $gatk CollectWgsMetrics \
       -I ${alignDeDup} \
       -O ${depthFile} \
       --READ_LENGTH 100 \
       -R ${refSeq}