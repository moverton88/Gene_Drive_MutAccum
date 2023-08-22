#!/bin/bash
#PBS -A mioverto
#PBS -l nodes=1 
#PBS -l walltime=7:00:00

while getopts 'a:R:s:' OPTION; do
   case "$OPTION" in
   v)
      variantFile="$OPTARG"
      echo "Variant file $OPTARG"
      ;;

   a)
      alignFile="$OPTARG"
      echo "Alignment file $OPTARG"
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

# variantFile=$sample
# $alignFile
# refSeq=$P1refSeq

refPrefix=${refSeq/.fna/}
refFile=$(basename "${refSeq}")
refName=${refFile%%_*}

if [ ! -f "${refPrefix}.fasta.img" ]; then

    if [ ! -f "${refPrefix}.fasta" ]; then
    echo "creating FASTA reference file"
    cp ${refSeq} ${refPrefix}.fasta
    fi
    echo "creating BWA index image file"
    java  -Xmx16G -jar $gatk BwaMemIndexImageCreator \
        -I ${refPrefix}.fasta \
        -O  ${refPrefix}.fasta.img

fi

# sampleName=$(basename "${variantFile}" .dedup.bam)
variantOut=${variantFile/.g.vcf/.faa.g.vcf}

java  -Xmx32G -jar $gatk FilterAlignmentArtifacts \
    -R ${refSeq} \
    -I ${alignFile} \
    --bwa-mem-index-image ${refPrefix}.fasta.img \
    --smith-waterman JAVA \
    -V ${variantFile} \
    -O ${variantOut}
