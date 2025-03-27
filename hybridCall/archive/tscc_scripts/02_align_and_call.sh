#!/bin/bash
#PBS -A mioverto
#PBS -l nodes=1 
#PBS -l walltime=24:00:00

while getopts 'p:r:v:t:R:s:' OPTION; do
   case "$OPTION" in
    p)
        pName="$OPTARG"
        echo "Parent name $OPTARG"
        ;;

    r)
        reads1="$OPTARG"
        echo "Reads file 1 $OPTARG"
        ;;

    v)
        gVCFdir="${OPTARG}"
        echo "VCF directory is $OPTARG"
        ;;
   
    t)
        VCFtype="${OPTARG}"
        echo "VCF type is $OPTARG"
        ;;

    R)
        refSeq="${OPTARG}"
        echo "Reference sequence $OPTARG"
        ;;

    s)
        opSys="${OPTARG}"
        echo "Operating system $OPTARG"
        ;;

    ?)
        echo "script usage: $(basename \$0) [-p parentName] [-r readsfile1] [-v VCFdirectory] [-t VCFtype] [-R referenceSequence]  [-s operatingSystem]"
    esac
done

# End Section #################################################################
###############################################################################
#                                                                             #
#                                                                             #
###############################################################################
###############################################################################
# script setup ################################################################

# if script arguemnets are not provided, exit
if [ $# -eq 0 ]; then
    >&2 echo "No arguments provided"
    exit 1
fi

# if running on TSCC, load modules, otherwise, assign paths to run applications
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

# End Section #################################################################
###############################################################################
#                                                                             #
#                                                                             #
###############################################################################
###############################################################################
# Set required variables ######################################################

# get reference sequence prefix
refPrefix=${refSeq/.fasta/}

# read 2 file path construction
reads2=${reads1/_R1/_R2}

# extract sample name from read file name
sampleName=$(basename "${reads1%%_R1*}")

# assign alignment paths and files
alignDirOut=${alignDir}/${pName}
alignRaw=${alignDirOut}/${sampleName}.bam
alignDeDup=${alignRaw/.bam/.dedup.bam}
deDupMetrics=${alignMetricsDir}/${pName}/${sampleName}_dupsMetric.txt
depthFile=${deDupMetrics/_dupsMetric.txt/_depth.txt}
tmpDir=${alignDirOut}/${sampleName}_tmp

# End Section #################################################################
###############################################################################
#                                                                             #
#                                                                             #
###############################################################################
###############################################################################
# If not present, create required directories #################################

if [ ! -d ${alignDirOut}/ ]; then
    echo "Create alignment output dir"
    mkdir ${alignDirOut}/
fi

if [ ! -d ${tmpDir} ]; then
    echo "Create alignment temp dir"
    mkdir ${alignDirOut}/${sampleName}_tmp
fi

if [ ! -d ${alignMetricsDir}/${pName}/ ]; then
    echo "Create metrics dir"
    mkdir ${alignMetricsDir}/${pName}/
fi

# End Section #################################################################
###############################################################################
#                                                                             #
#                                                                             #
###############################################################################
# Indexed and pre-process reference sequence if not already done ##############

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

# End Section #################################################################
###############################################################################
#                                                                             #
#                                                                             #
###############################################################################
# Align reads with bwa ########################################################

# ${sampleName} used for naming alignments. Output is a sorted bam alignment.
# bwa mem -v [verbose level] 
   # align reads with default bwa settings
# samtools view -b
   # convert alignment from .sam to .bam
# samtools sort -o [output file] -T [temp files] -o [outputFile]
   # sort the .bam file according to position
# samtools index [inputFile]
   # Create bam index file
# bwa mem -v 3 $refPrefix $reads1 $reads2
bwa mem -v 3 $refPrefix $reads1 $reads2 | samtools view -b - | \
samtools sort -O bam -T ${tmpDir} -o ${alignRaw} - 
samtools index ${alignRaw}

rm -dr ${alignDirOut}/${sampleName}_tmp

# End Section #################################################################
###############################################################################
#                                                                             #
#                                                                             #
###############################################################################
# Standardize read group labels ###############################################

## GATK requres a different output file, cannot overwrite existing bam file
## --RGPL [platform]
## --RGPU [read group platform unit]
## --RGLB [read group library]
## --RGSM [read group sample name]

java -jar $gatk AddOrReplaceReadGroups \
   -I ${alignRaw} \
   -O ${alignRaw/.bam/.rg.bam} \
   --RGPL ILLUMINA \
   --RGPU ${proj}_${sampleName} \
   --RGLB ${proj}_${pName} \
   --RGSM ${sampleName}

# index bam alignment
samtools index ${alignRaw/.bam/.rg.bam}

# remove initial bam files
rm ${alignRaw}
rm ${alignRaw/.bam/.bam.bai}

# End Section #################################################################
###############################################################################
#                                                                             #
#                                                                             #
###############################################################################
###############################################################################
# Remove duplicate reads and create metrics file
java -jar $gatk MarkDuplicates \
   -I ${alignRaw/.bam/.rg.bam} \
   -O ${alignDeDup} \
   -M ${deDupMetrics} \
   --TMP_DIR ${alignDir}/${pName}/${sampleName}_tmp

# remove raw alignment files
rm ${alignRaw/.bam/.rg.bam}
rm ${alignRaw/.bam/.rg.bam.bai}

# index bam alignment
samtools index ${alignDeDup}

# Generate depth metrics file
java -jar -Xmx16g $gatk CollectWgsMetrics \
       -I ${alignDeDup} \
       -O ${depthFile} \
       --READ_LENGTH 100 \
       -R ${refSeq}

# End Section #################################################################
###############################################################################
#                                                                             #
#                                                                             #
###############################################################################
###############################################################################
# Call variants with GATK

# Depending on VCF type, emits either a normal VCF for immediate use or a gVCF,
# which is combined with all gVCFs of a cohort and called using GenotypeGVCF

if [ ! -d "${gVCFdir}" ]; then
    echo "Create gVCF variants dir"
    mkdir ${gVCFdir}
fi

# Call variants with HaplotypeCaller. -ERC generates a gVCF file

if [[ ${VCFtype} == "g" || ${VCFtype} == "gVCF" || ${VCFtype} == "gvcf" ]]; then

gVCFout=${gVCFdir}/${sampleName}.g.vcf

java  -Xmx16G -jar $gatk HaplotypeCaller \
    -R ${refSeq} \
    -I ${alignDeDup} \
    -O ${gVCFout} \
    --read-filter MappingQualityReadFilter \
    --minimum-mapping-quality 30 \
    --max-reads-per-alignment-start 100 \
    --max-assembly-region-size 500 \
    -ERC GVCF 

elif [[ ${VCFtype} == "v" || ${VCFtype} == "VCF" || ${VCFtype} == "vcf" ]]; then

gVCFout=${gVCFdir}/${sampleName}.vcf

java  -Xmx16G -jar $gatk HaplotypeCaller \
    -R ${refSeq} \
    -I ${alignDeDup} \
    -O ${gVCFout} \
    --read-filter MappingQualityReadFilter \
    --minimum-mapping-quality 30 \
    --max-reads-per-alignment-start 100 \
    --max-assembly-region-size 500

else

   echo "VCF output type not recognized"
fi
