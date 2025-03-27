#!/bin/bash

###############################################################################
# This script assigns all of the variables that might be needed for any of the
# pipeline processes. As such, after initializing the project, you must place
# The required files in the appropriate directories created during 
# initialization. This includes:
#   - All external applications -> [homeDir]/bin/
#   - The type reference sequence -> ./references/
#   - The parent accession list file -> ./metadata/
#   - The clone accession list file  -> ./metadata/
#   - The parental rename file -> ./metadata/
#   - The clone rename file -> ./metadata/

###############################################################################

projDir=$1
echo "Project directory is set to $projDir"
typeRef=$2
echo "The name of the type reference sequence is $typeRef"
P1=$3
echo "Parent 1 name is $P1"
P2=$4
echo "Parent 2 name is $P2"

###############################################################################
# Read the sysparams.txt file for required path locations
export opSys=`sed -n 1p ${projDir}/sysparams.txt`
export proj=`sed -n 2p ${projDir}/sysparams.txt`
export appDir=`sed -n 4p ${projDir}/sysparams.txt`
export codeDir=`sed -n 5p ${projDir}/sysparams.txt`

export sysVars=${projDir}/sysvars.sh
source $sysVars

# Trimmomatic
export bwa=${appDir}/bwa-mem2-2.2.1_x64-linux/bwa-mem2
export samtools=${appDir}/samtools-1.19.2/samtools
export bcftools=${appDir}/bcftools-1.19/bcftools
export bedtools=${appDir}/bedtools
export bgzip=${appDir}/htslib-1.19.1/bgzip
# export java=${appDir}/jdk17/bin/java
export java=${appDir}/jdk-8u402-b06-linux-x64/bin/java
export gatk=${appDir}/gatk-4.3.0.0/gatk-package-4.3.0.0-local.jar
export nucmer=${appDir}/mummer-4.0.0rc1/nucmer
export delta=${appDir}/mummer-4.0.0rc1/delta-filter
export coords=${appDir}/mummer-4.0.0rc1/show-coords
export delta2vcf=${appDir}/mummer-4.0.0rc1/delta2vcf
export fastqc=${appDir}/FastQC/fastqc

## set all variables used in pipeline. Need to ensure alternatives to TSCC modules
export refDir=${projDir}/references
export refSeqDir=${refDir}/final
export typeRefDir=`find ${refSeqDir} -type d -name ${typeRef}`
export typeRefSeq=`find ${typeRefDir} -type f -name *${typeRef}*.fasta`

export constDir=${refDir}/construction
export pReadsDir=${constDir}/reads
export pReads1=${pReadsDir}/${P1}_R1.fastq
export pReads2=${pReadsDir}/${P1}_R2.fastq
export pAlignDir=${constDir}/alignments
export pDeDupMetrics=${constDir}/metrics
export pVarDir=${constDir}/variants


if [[ ! -d ${refSeqDir}/${P1}/ ]]; then
    # echo "P1 reference dir missing"
    mkdir ${refSeqDir}/${P1}/
fi
if [[ ! -d ${refSeqDir}/${P2}/ ]]; then
    # echo "P1 reference dir missing"
    mkdir ${refSeqDir}/${P2}/
fi

export P1refSeq=${refSeqDir}/${P1}/${P1}_refseq.fasta
export P2refSeq=${refSeqDir}/${P2}/${P2}_refseq.fasta

export readsDir=${projDir}/reads
export alignDir=${projDir}/alignments
export variantsDir=${projDir}/variants
export indvVarDir=${variantsDir}/individuals
export lineVarDir=${variantsDir}/lines
export finalVarDir=${variantsDir}/final
export readMetricsDir=${projDir}/metrics/reads
export alignMetricsDir=${projDir}/metrics/alignments
export metaDir=${projDir}/metadata

export accFile=`find $metaDir -type f -name *clone*ccession*`
export pAccFile=`find $metaDir -type f -name *arent*ccession*`
export renameFile=`find $metaDir -type f -name *clone*rename.*`
export pRenameFile=`find $metaDir -type f -name *parent*rename.*`
# export trimFile=`find $metaDir -type f -name *[Tt]rim*`
