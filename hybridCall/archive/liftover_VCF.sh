#!/bin/bash


while getopts 'v:c:R:' OPTION; do
   case "$OPTION" in
   v)
      vcfIn="$OPTARG"
      echo "Reads file 1 $OPTARG"
      ;;

   c)
      chain="${OPTARG}"
      echo "Operating system $OPTARG"
      ;;

    R)
      typeRefSeq="$OPTARG"
      echo "Reads file 1 $OPTARG"
      ;;

   ?)
        echo "script usage: $(basename \$0) [-R1 readsfile1] [-s setvarablesfile]"
    esac
done

# P1=BTI
# vcfIn=${finalVarDir}/${P1}/all_${P1}.final.vcf
# chain=${refDir}/final/${P1}/${P1}xS288C.chain
# typeRefSeq=${refDir}/final/S288C/S288C_R64_refseq.fna 

java -jar $gatk LiftoverVcf \
    -I ${vcfIn} \
    -O ${vcfIn/.vcf/.lift.vcf} \
    -C $chain \
    --REJECT ${vcfIn/.vcf/.reject.vcf} \
    -R $typeRefSeq

