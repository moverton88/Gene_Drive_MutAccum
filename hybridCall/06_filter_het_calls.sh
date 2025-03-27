#!/bin/bash
#PBS -A mioverto
#PBS -l nodes=1 
#PBS -l walltime=24:00:00


while getopts 's:m:t:r:S:' OPTION; do
   case "$OPTION" in
    s)
        SNPs_input="$OPTARG"
        echo "Reconciled SNPs table $OPTARG"
        ;;

    m)
        model_input="$OPTARG"
        echo "Model parameters list $OPTARG"
        ;;

    t)
        model_type="${OPTARG}"
        echo "Type of het model to use ('v' = modeled variance, 'd' or 'n' = median or mean BB alpha): $OPTARG"
        ;;

    r)
        FE_rate="${OPTARG}"
        echo "Rate of falsely excluded true het calls $OPTARG"
        ;;

    S)
        SNPs_output="${OPTARG}"
        echo "SNPs table out $OPTARG"
        ;;

    ?)
        echo "script usage: $(basename \$0) [-s SNPsIn] [-m modelInput] [-t modelTypeVarianceorAlpha] [-r FalseExsclusionRate] [-S SNPsOut]"
    esac
done

# if script arguemnets are not provided, exit
if [ $# -eq 0 ]; then
    >&2 echo "No arguments provided"
    exit 1
fi

# samples="all"
# SNPs_input=${variantsDir}/reconciled/${samples}_${P1}x${P2}_clean_SNPs.RData
# model_input=${variantsDir}/reconciled/${samples}_BYxRM_hetTrain.RData
# model_type="v"
# FE_rate=0.05
# SNPs_output=${variantsDir}/reconciled/${samples}_BYxRM_hetFilter_SNPs.RData

# End Section #################################################################
###############################################################################

export filter_hets=/home/mioverto/code/hybridCall/06b_filter_het_calls.R

module load R
Rscript $filter_hets -s $SNPs_input -m $model_input -t $model_type -r $FE_rate -S $SNPs_output
