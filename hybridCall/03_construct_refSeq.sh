#!/bin/bash
#PBS -A mioverto
#PBS -l nodes=1 
#PBS -l walltime=24:00:00

###############################################################################
# 03_construct_refSeq.sh
# Takes in a set of paired read sequence files, aligns them to a Type
# reference, calls SNP and indel variants, and constructs a consensus sequence
# and a liftover file.
#*****************************************************************************#


while getopts 'R:s:' OPTION; do
   case "$OPTION" in
   R)
      pReads1="$OPTARG"
      echo "Reads file 1 $OPTARG"
      ;;

   s)
      opSys="${OPTARG}"
      echo "Operating system $OPTARG"
      ;;
   ?)
        echo "script usage: $(basename \$0) [-R1 readsfile1] [-s operatingSystem]"
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

###############################################################################
# Manually entered variables
#*****************************************************************************#
# pReads1=${pReadsDir}/${P1}*R1*fastq
# pReads1=`find ${pReadsDir} -type f -name "${P1}*R1*fastq"`
# pAlignment=${pAlignDir}/${pName}.bam
###############################################################################

###############################################################################
# Set required variables
#*****************************************************************************#

pReads2=${pReads1/_R1./_R2.}
# Get Parent name from filename
pFileName=$(basename "${pReads1}")
pName=${pFileName%%_*}

# Get name of Type strain from filename
typeRefPrefix=${typeRefSeq/.fasta/}
typefn=$(basename "${typeRefSeq}" .fasta)
typeName=${typefn%%_*}

# Alignment and variant path names
pAlignment=${pAlignDir}/${pName}.bam
pVariant=${pVarDir}/${pName}.vcf

# Path and file names for final parental reference sequence and liftover files
pReference=${refSeqDir}/${pName}/${pName}_refseq.fasta
pChain=${refSeqDir}/${pName}/${pName}x${typeName}.chain

###############################################################################


###############################################################################
# Create required files if absent
## BWA index file of Type reference sequence
## GATK sequence dictionary of Type reference sequence
## faidx index of Type reference sequence
#*****************************************************************************#

if [[ -f "${typeRefPrefix}.amb" ]]; then
   echo "bwa reference index files already exists"
else
   echo "bwa reference index files does not exist. Generating now..."   
   bwa index -p $typeRefPrefix $typeRefSeq 
fi

if [[ -f "${typeRefSeq}.fai" ]]; then
   echo "GATK reference index file already exists"
else
   echo "GATK reference index file does not exist. Generating now..."   
   java -jar $gatk CreateSequenceDictionary \
      -R $typeRefSeq \
      -O ${typeRefSeq/.fasta/.dict}

   samtools faidx $typeRefSeq
fi

###############################################################################

###############################################################################
# Align parental reference with type reference using bwa
# Convert sam to bam, sort by positions, and index
#*****************************************************************************#

bwa mem -v 1 $typeRefPrefix $pReads1 $pReads2 | samtools view -b - | \
samtools sort -O bam -o ${pAlignment} - 

samtools index ${pAlignment}

###############################################################################

###############################################################################
# bwa does not correctly assign read groups, so must assign manually
#*****************************************************************************#

java -jar $gatk AddOrReplaceReadGroups \
   -I ${pAlignment}  \
   -O ${pAlignment/.bam/.crct.bam} \
   -R $typeRefSeq \
   -LB lib1 \
   -PL ILLUMINA \
   -PU assembly \
   -SM ${pName}

###############################################################################

###############################################################################
# Remove duplicate reads and re-index
#*****************************************************************************#

java -jar $gatk MarkDuplicates \
   --INPUT ${pAlignment/.bam/.crct.bam} \
   --OUTPUT ${pAlignment/.bam/.dedup.bam} \
   --METRICS_FILE ${pDeDupMetrics}/${pName}_dedup.txt \
   --ASSUME_SORTED TRUE \
   --TMP_DIR ${pAlignDir}/${pName}_tmp

# rm ${pAlignment}
# rm ${pAlignment/.bam/.bam.bai}

samtools index ${pAlignment/.bam/.dedup.bam}

###############################################################################

###############################################################################
# Call all variants with haplotype caller
#*****************************************************************************#

java -Xmx16G -jar $gatk HaplotypeCaller \
   -R $typeRefSeq \
   -I ${pAlignment/.bam/.dedup.bam} \
   -O ${pVariant}

###############################################################################

###############################################################################
# Separate SNP and indel variants and filter each with GATK rec'd settings
#*****************************************************************************#

java -jar -Xmx16g $gatk SelectVariants  \
    -R $typeRefSeq \
    --variant ${pVariant} \
    --select-type-to-include INDEL \
    -O ${pVariant/.vcf/.indel.vcf}
    
java -jar -Xmx16g $gatk SelectVariants  \
    -R $typeRefSeq \
    --variant ${pVariant} \
    --select-type-to-include SNP \
    -O ${pVariant/.vcf/.snps.vcf}

rm ${pVariant}

# SNP filter
java -jar $gatk VariantFiltration \
    -R $typeRefSeq \
    -V ${pVariant/.vcf/.snps.vcf} \
    -O ${pVariant/.vcf/.snps.m.vcf} \
    --filter-name "QD2" --filter-expression "QD < 2.0" \
    --filter-name "FS60" --filter-expression "FS > 60.0" \
    --filter-name "MQ20" --filter-expression "MQ < 20.0" \
    --filter-name "QUAL30" --filter-expression "QUAL < 30.0" \
    --filter-name "SOR10" --filter-expression "SOR > 10.0" \
    --filter-name "DP6" --filter-expression "DP < 6"

java -jar $gatk SelectVariants \
    -V ${pVariant/.vcf/.snps.m.vcf} \
    -O ${pVariant/.vcf/.snps.f.vcf} \
    --exclude-filtered true

# Indel filter
java -jar $gatk VariantFiltration \
    -R $typeRefSeq \
    -V ${pVariant/.vcf/.indel.vcf} \
    -O ${pVariant/.vcf/.indel.m.vcf} \
    --filter-name "QD2" --filter-expression "QD < 2.0" \
    --filter-name "FS200" --filter-expression "FS > 200.0" \
    --filter-name "QUAL30" --filter-expression "QUAL < 30.0" \
    --filter-name "DP6" --filter-expression "DP < 6"

java -jar $gatk SelectVariants \
    -V ${pVariant/.vcf/.indel.m.vcf} \
    -O ${pVariant/.vcf/.indel.f.vcf} \
    --exclude-filtered true 

###############################################################################

###############################################################################
# Merge filtered vcfs and filter each call by GQ likelihood score >= 30 and 
# GT = 1/1
#*****************************************************************************#

java -jar -Xmx16g $gatk MergeVcfs  \
   -I ${pVariant/.vcf/.snps.f.vcf} \
   -I ${pVariant/.vcf/.indel.f.vcf}  \
   -O ${pVariant}

rm ${pVarDir}/${pName}*indel*
rm ${pVarDir}/${pName}*snps*

java -jar $gatk FilterVcf \
   -I ${pVariant} \
   --MIN_GQ 30 \
   -O ${pVariant/.vcf/.inter.vcf}

bcftools view -i 'GT="1/1"' -o ${pVariant} ${pVariant/.vcf/.inter.vcf}

###############################################################################

###############################################################################
# Create consensus sequence and liftover files from vcf file and 
# Type reference sequence
#*****************************************************************************#

bgzip -c ${pVariant} > ${pVariant/.vcf/.vcf.gz}
bcftools index ${pVariant/.vcf/.vcf.gz}
bcftools consensus -H A -f ${typeRefSeq} -c $pChain -o $pReference ${pVariant/.vcf/.vcf.gz} 

rm ${pVarDir}/*inter*

###############################################################################
