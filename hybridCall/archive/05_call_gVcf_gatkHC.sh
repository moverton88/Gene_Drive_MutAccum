#!/bin/bash
#PBS -A mioverto
#PBS -l nodes=1 
#PBS -l walltime=7:00:00

while getopts 'a:v:R:t:s:' OPTION; do
   case "$OPTION" in
   a)
      alignFile="$OPTARG"
      echo "Alignment file $OPTARG"
      ;;

   v)
      gVCFdir="${OPTARG}"
      echo "VCF directory is $OPTARG"
      ;;

   R)
      refSeq="${OPTARG}"
      echo "Reference sequence $OPTARG"
      ;;
   
   t)
      VCFtype="${OPTARG}"
      echo "VCF type is $OPTARG"
      ;;

   s)
      opSys="${OPTARG}"
      echo "Operating system $OPTARG"
      ;;

   ?)
        echo "script usage: $(basename \$0) [-a alignmentFile] [-v gVCFdirectory] [-R referenceSequence] [-t VCFtype] [-s operatingSystem]"
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

# Start with alignment. Assumes renamed, trimmed fastq files and un-indexed reference fasta file

# alignFile=${alignDir}/${P1}/${line}_${P1}.dedup.bam
# refSeq=$P1refSeq
# refSeq=${P1refSeq/.fna/_mum.fna}
# alignFile=${alignDir}/${refName}/${sampleName}x${refName}.dedup.bam
# sampleName=$(basename "${alignFile}" .dedup.bam)
# gVCFdir=${variantsDir}/individuals/${P1}

refPrefix=${refSeq/.fasta/}
refFile=$(basename "${refSeq}")
refName=${refFile%%_*}

if [ ! -d "${gVCFdir}" ]; then
    echo "Create gVCF variants dir"
    mkdir ${gVCFdir}
fi

# Call variants with HaplotypeCaller. -ERC generates a gVCF file

if [[ ${VCFtype} == "g" || ${VCFtype} == "gVCF" || ${VCFtype} == "gvcf" ]]; then

gVCFout=${gVCFdir}/${sampleName}.g.vcf

java  -Xmx16G -jar $gatk HaplotypeCaller \
    -R ${refSeq} \
    -I ${alignFile} \
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
    -I ${alignFile} \
    -O ${gVCFout} \
    --read-filter MappingQualityReadFilter \
    --minimum-mapping-quality 30 \
    --max-reads-per-alignment-start 100 \
    --max-assembly-region-size 500

else

   echo "VCF output type not recognized"
fi
