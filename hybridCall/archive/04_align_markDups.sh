#!/bin/bash
#PBS -A mioverto
#PBS -l nodes=1 
#PBS -l walltime=24:00:00

while getopts 'r:a:R:s:' OPTION; do
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

   s)
      opSys="$OPTARG"
      echo "Operating system $OPTARG"
      ;;

   ?)
        echo "script usage: $(basename \$0) [-r readsfile1] [-a alignmentDirOut] [-R referenceSequence]  [-s operatingSystem]"
    esac
done


if [[ ${opSys} == "tscc" ]]; then
   echo "loading packages from modules"
      module load bwa
      module load samtools    
      module load bcftools
      module load bedtools
   else
   echo "loading packages from files"
      bwa=${bwaApp}
      samtools=${samtoolsApp}
      bcftools=${bcftoolsApp}
      bedtools=${bedtoolsApp}
fi

# refSeq=$P1refSeq
# refSeq=${P1refSeq/.fasta/_mum.fasta}
refPrefix=${refSeq/.fasta/}
refFile=$(basename "${refSeq}")
refName=${refFile%%_*}

# reads1=${readsFile}
reads2=${reads1/_R1/_R2}
# reads2=${reads1/_R1P/_R2P}
readsName=$(basename "${reads1%%_R1*}")

if [ -z "$sampleName" ]; then  
   sampleName=${readsName}
fi

# alignDirOut=${alignDir}/${P1}
alignRaw=${alignDirOut}/${sampleName}.bam
alignDeDup=${alignRaw/.bam/.dedup.bam}
deDupMetrics=${alignMetricsDir}/${refName}/${sampleName}_dupsMetric.txt
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

if [ ! -d ${alignMetricsDir}/${refName}/ ]; then
    echo "Create metrics dir"
    mkdir ${alignMetricsDir}/${refName}/
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

if [[ ! -f "${refPrefix}.dict" || ! -f "${refPrefix}.fasta.fai" ]]; then
   
   echo "GATK reference index file does not exist. Generating now..."   
   java -jar $gatk CreateSequenceDictionary \
      -R ${refSeq} \
      -O ${refPrefix}.dict

   samtools faidx $refSeq

else
   echo "GATK reference index file already exists"
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

bwa mem -v 3 $refPrefix $reads1 $reads2 | samtools view -b - | \
samtools sort -O bam -T ${tmpDir} -o ${alignRaw} - 
samtools index ${alignRaw}

rm -dr ${alignDirOut}/${sampleName}_tmp

### Standardize read group labels
# GATK requres a different output file, cannot overwrite existing bam file
# --RGPL [platform]
# --RGPU [read group platform unit]
# --RGLB [read group library]
# --RGSM [read group sample name]


java -jar $gatk AddOrReplaceReadGroups \
   -I ${alignRaw} \
   -O ${alignRaw/.bam/.rg.bam} \
   --RGPL ILLUMINA \
   --RGPU ${proj}_${sampleName} \
   --RGLB ${proj}_${refName} \
   --RGSM ${readsName}

# index bam alignment
samtools index ${alignRaw/.bam/.rg.bam}

# remove initial bam files
rm ${alignRaw}
rm ${alignRaw/.bam/.bam.bai}


# Remove duplicate reads and create metrics file
java -jar $gatk MarkDuplicates \
   -I ${alignRaw/.bam/.rg.bam} \
   -O ${alignDeDup} \
   -M ${deDupMetrics} \
   --TMP_DIR ${alignDir}/${refName}/${sampleName}_tmp

# remove raw alignment files
rm ${alignRaw/.bam/.rg.bam}
rm ${alignRaw/.bam/.rg.bam.bai}

# index bam alignment
samtools index ${alignDeDup}

# rm -r ${bamDir}/temp_dir/${index}
# rm ${alignRaw}
# rm ${alignRaw/.bam/.bam.bai}

# Generate depth metrics file
java -jar -Xmx16g $gatk CollectWgsMetrics \
       -I ${alignDeDup} \
       -O ${depthFile} \
       --READ_LENGTH 100 \
       -R ${refSeq}
