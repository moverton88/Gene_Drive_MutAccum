#!/bin/bash
#PBS -A mioverto
#PBS -l nodes=1 
#PBS -l walltime=5:00:00

# Combine VCF files of individual clones into one master VCF
while getopts 'v:l:f:p:R:' OPTION; do
   case "$OPTION" in
    v)
        variantDir="$OPTARG"
        echo "Variant directory $OPTARG"
    ;;

    l)
        lineage="$OPTARG"
        echo "Lineage $OPTARG"
    ;;

    f)
        founders="$OPTARG"
        if [[ -f ${founders} ]]; then
            echo "The founders of this lineage are \n$(cat $founders)\n"
        elif [[ -z ${founders} || ${founders} == "" ]]; then
            echo "No founders provided"
            founders="none"
        else
            echo "The founder of this lineage is ${founders}"
        fi
    ;;

   p)
        pName="${OPTARG}"
        echo "Parent name is $OPTARG"
    ;;

   R)
        refSeq="${OPTARG}"
        echo "reference sequence $OPTARG"
    ;;

   ?)
        echo "script usage: $(basename \$0) [-v variantDirectory] [-l lineage] \
         [-f founderNameOrFileOfNames] [-p parentName] [-R referenceSequence]"
    esac
done

# founders=none
# variantDir=${variantsDir}/individuals/${P1}
# lineage=${line}
# refSeq=$P1refSeq
refFile=$(basename "${refSeq}")


if [ ! -d "${lineVarDir}/${pName}" ]; then
    echo "Creating lineage variants dir"
    mkdir ${lineVarDir}/${pName}
fi

if [ ! -d "${finalVarDir}/${pName}" ]; then
    echo "Creating final variants dir"
    mkdir ${finalVarDir}/${pName}
fi


export gVCFlist=${lineVarDir}/${pName}/${lineage}_gVCFlist.list

if [[ ${lineage} == all || ${lineage} == All ]] ; then
    echo "making list of all line samples"
    dir ${variantDir}/*${pName}*.g.vcf > ${gVCFlist}

else
    echo "making list of samples from ${lineage}"
    find ${variantDir} -name "${lineage}*.g.vcf" > ${gVCFlist}

fi

lineVCF=${lineVarDir}/${pName}/${lineage}_${pName}.line.vcf
finalVCF=${finalVarDir}/${pName}/${lineage}_${pName}.final.vcf

# export gVCFlist=${variantDir}/anc_gVCFlist.list
# dir ${variantDir}/*00*.g.vcf >> ${gVCFlist}
# less ${gVCFlist}
# less ${founderList}


if [[ ${founders} == "none" ]]; then
    echo "no founder IDs provided, genotyping end-point samples only"
    java -jar -Xmx16g $gatk CombineGVCFs  \
        -R ${refSeq} \
        -V ${gVCFlist} \
        -O ${lineVCF}

    java -jar -Xmx16g $gatk GenotypeGVCFs \
        -R ${refSeq} \
        -V ${lineVCF} \
        -O ${finalVCF}

else

    echo "founder IDs provided, genotyping all samples"
    java -jar -Xmx16g $gatk CombineGVCFs  \
        -R ${refSeq} \
        -V $gVCFlist \
        --founder-id ${founders} \
        -O ${lineVCF}

    java -jar -Xmx16g $gatk GenotypeGVCFs \
        -R ${refSeq} \
        -V ${lineVCF} \
        --founder-id ${founders} \
        -O ${finalVCF}
        
fi
