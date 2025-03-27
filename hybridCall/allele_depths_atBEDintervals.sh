#!/bin/bash
#PBS -A mioverto
#PBS -l nodes=1 
#PBS -l walltime=5:00:00


while getopts 'p:a:b:o:R:' OPTION; do
   case "$OPTION" in
    p)
        pName="$OPTARG"
        echo "Parent name $OPTARG"
        ;;

    a)
        align_in="$OPTARG"
        echo "Bam alignment file $OPTARG"
        ;;

    b)
        bed_in="$OPTARG"
        echo "BED intervals file $OPTARG"
        ;;

    o)
        tsv_out="${OPTARG}"
        echo "Table file out $OPTARG"
        ;;

    R)
        refSeq="${OPTARG}"
        echo "Reference sequence $OPTARG"
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

# align_in=${alignDir}/${P1}/N_A00.dedup.bam
# bed_in=${metaDir}/hom_sites.bed
# tsv_out=${alignDir}/${P1}/N_A00_HC.tsv
# refSeq=${refDir}/final/${P1}/${P1}_refseq.fasta


java  -Xmx16G -jar $gatk HaplotypeCaller \
    -R ${refSeq} \
    -I ${align_in} \
    -O ${tsv_out/.tsv/.vcf} \
    -L ${bed_in} \
    --read-filter MappingQualityReadFilter \
    -ERC BP_RESOLUTION \
    --output-mode EMIT_ALL_CONFIDENT_SITES

java  -Xmx16G -jar $gatk VariantsToTable \
     -V ${tsv_out/.tsv/.vcf} \
     -F CHROM -F POS -GF GT -GF AD -GF GQ \
     -O ${tsv_out}
