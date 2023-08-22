#!/bin/bash
#PBS -A mioverto
#PBS -l nodes=1 
#PBS -l walltime=24:00:00

###############################################################################
# 03_construct_refSeq_contigs.sh
# Takes in a set of scaffold or contig-level sequences, aligns them to a Type
# reference, calls SNP and indel variants, and constructs a consensus sequence
# and a liftover file.
#*****************************************************************************#

while getopts 'R:s:' OPTION; do
   case "$OPTION" in
   R)
      pReads="$OPTARG"
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
# Requires the Mummer program from 
# https://mummer4.github.io/index.html
#*****************************************************************************#
nucmer=${appDir}/mummer-4.0.0rc1/nucmer
delta=${appDir}/mummer-4.0.0rc1/delta-filter
# coords=${appDir}/mummer-4.0.0rc1/show-coords
delta2vcf=${appDir}/mummer-4.0.0rc1/delta2vcf

###############################################################################
# Manually entered variables
#*****************************************************************************#
# pReads=${pReadsDir}/${P1}_contigs.fasta
# pReads=`find ${pReadsDir} -type f -name "${P1}*contigs.fasta"`
# pAlignment=${pAlignDir}/${pName}.bam
###############################################################################

###############################################################################
# Set required variables
#*****************************************************************************#

# Get Parent name from filename
pFileName=$(basename "${pReads}")
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
# BWA index file of Type reference sequence
# GATK sequence dictionary of Type reference sequence
# faidx index of Type reference sequence
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
# Map scaffold sequences to Type sequence using nucmer from the Mummer program
# Call variants with Mummer delta function
# Correct GT field in vcf file
#*****************************************************************************#

cd ${pAlignDir}/${pName} # nucmer only outputs to stdout so must change to output dir
$nucmer -p ${pName}_mum_ref ${typeRefSeq} ${pReads} # -p prefix of output
$delta ${pName}_mum_ref.delta -1 | $delta2vcf > ${pName}_mum_ref.vcf

# nucmer outputs the GT field as an integer rather than 0/0 ect format, so have to change
sed 's/0:1:/1\/1:1:/g' ${pName}_mum.vcf | \
sed 's/GT,Number=1,Type=Integer/GT,Number=1,Type=String/g' > ${pVariant}

###############################################################################

###############################################################################
# To create a consensus sequence, bcftools requires the vcf file to be 
# normalized. However, if there are overlapping indel and SNP variants, 
# the consensus call is not clear. 
#*****************************************************************************#

firstNorm=${pName}/.nrm.bcf.gz
SNPvariants=${pName}.snps.bcf.gz
SNPconsSeq=${pName}_SNPs.fasta
indelVariants=${pName}.indels.bcf.gz
indelNormVariants=${pName}.nrmIndels.bcf.gz

# First pass normalization of all variants
bcftools norm -m +any -f $typeRefSeq -Ob ${pVariant} > ${firstNorm}
bcftools index ${pName}.nrm.bcf.gz

# Split off SNP variants and create SNP consensus sequence
bcftools view -v snps -Ob ${firstNorm} -o ${SNPvariants}
bcftools index ${SNPvariants}

# Create SNP consensus sequence from 1/1 genotype SNPs
bcftools consensus -H A -f ${typeRefSeq} ${SNPvariants} > ${SNPconsSeq}

# Split off indel variants and normalize against SNP consensus sequence
bcftools view -v indels -Ob ${firstNorm} -o ${indelVariants}
bcftools norm -m +any -c s -f ${SNPconsSeq} -Ob ${indelVariants} > ${indelNormVariants}
bcftools index ${indelNormVariants}

# Create final consensus sequence by appending SNP consensus with 1/1 indel variants
bcftools consensus -H A -f ${SNPconsSeq} -c ${pChain} ${indelNormVariants} > ${pReference}

###############################################################################