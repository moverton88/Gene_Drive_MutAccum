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
      refSeq="${OPTARG}"
      echo "Reference sequence $OPTARG"
      ;;

   s)
      setVars="${OPTARG}"
      echo "Operating system $OPTARG"
      ;;

   ?)
        echo "script usage: $(basename \$0) [-R1 readsfile1] [-s setvarablesfile]"
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


# Reference must be indexed by bowtie2 before read alignment
# if [ ! -f ${refPrefix}.1.bt2 ]; then
#     echo "Reference sequence needs to be indexed (bowtie2-build)"
#     bowtie2-build ${refSeq} ${refPrefix}
# fi

if [[ ! -f "${refPrefix}.amb" || ! -f "${refPrefix}.bwt" || ! -f "${refPrefix}.sa" ]]; then
   echo "bwa reference index files does not exist. Generating now..."   
   bwa index -p $refPrefix $refSeq

else
   echo "bwa reference index files already exists"

fi

if [[ ! -f "${refSeq}.img" ]]; then
   echo "bwa index image file does not exist. Generating now..."   

   # cp ${refSeq} ${refPrefix}.fasta

   java -jar $gatk BwaMemIndexImageCreator \
     -I ${refPrefix}.fasta \
     -O ${refSeq}.img
   
   # rm ${refPrefix}.fasta

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

# Align paired and unpaired reads reads with bowtie2. 
# ${sampleName} used for naming alignments. Output is sorted bam alignment
# -I, -X flag give the minimum, maximum valid gap between paired reads
# samtools view -h [include header] -b [bam output] -u [uncompressed .bam] | sort -m [memory allocated] -o [output file] -T [temp files]

# bowtie2 --rg-id ${sampleName} --rg SM:${sampleName} --sensitive-local -I 100 -X 1000 \
# -x ${refPrefix} -1 ${R1PFILE} -2 ${R2PFILE} -U "${R1UFILE},${R2UFILE}" \
# | samtools view -h -b -u | samtools sort -m 10000000 -o ${bamRaw} -O bam -T ${bamDir}/temp_dir/${sampleName}

bwa mem -v 3 $refPrefix $reads1 $reads2 | samtools view -b - | \
samtools sort -O bam -T ${tmpDir} -o ${alignRaw} - 
samtools index ${alignRaw}

rm -dr ${alignDirOut}/${sampleName}_tmp

### Standardize read group labels

java -jar $gatk AddOrReplaceReadGroups \
   -I ${alignRaw} \
   -O ${alignRaw/.bam/.rg.bam} \
   --RGPL ILLUMINA \
   --RGPU ${proj}_${sampleName} \
   --RGLB ${proj}_${refName} \
   --RGSM ${readsName}

# index bam alignment
samtools index ${alignRaw/.bam/.rg.bam}

rm ${alignRaw}
rm ${alignRaw/.bam/.bam.bai}

# other bowtie2 options: -N 1 -L 20 --mp 3,2 --rdg 3,2 --rfg 3,2 --local
# Remove duplicate reads
java -jar $gatk MarkDuplicates \
   -I ${alignRaw/.bam/.rg.bam} \
   -O ${alignDeDup} \
   -M ${deDupMetrics} \
   --TMP_DIR ${alignDir}/${refName}/${sampleName}_tmp

rm ${alignRaw/.bam/.rg.bam}
rm ${alignRaw/.bam/.rg.bam.bai}

# index bam alignment
samtools index ${alignDeDup}

# rm -r ${bamDir}/temp_dir/${index}
# rm ${alignRaw}
# rm ${alignRaw/.bam/.bam.bai}

java -jar -Xmx16g $gatk CollectWgsMetrics \
       -I ${alignDeDup} \
       -O ${depthFile} \
       --READ_LENGTH 100 \
       -R ${refSeq} 

