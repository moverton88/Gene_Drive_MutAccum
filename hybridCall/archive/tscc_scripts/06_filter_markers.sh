#!/bin/bash
#PBS -A mioverto
#PBS -l nodes=1 
#PBS -l walltime=24:00:00


while getopts 's:t:c:h:r:S:' OPTION; do
   case "$OPTION" in
    s)
        SNPs_merge="$OPTARG"
        echo "Reconciled SNPs table $OPTARG"
        ;;

    t)
        train_ID="$OPTARG"
        echo "IDs for model training data $OPTARG"
        ;;

    f)
        fractionTrainIDs="${OPTARG}"
        echo "Fraction of training calls for site to be known Het: $OPTARG"
        ;;
   
    F)
        fractionNonTrainIDs="${OPTARG}"
        echo "Fraction of non-training calls for site to be known Het: $OPTARG"
        ;;

    r)
        repeats_bed="${OPTARG}"
        echo "Repeats BED file $OPTARG"
        ;;

    S)
        SNPs_output="${OPTARG}"
        echo "SNPs table out $OPTARG"
        ;;

    ?)
        echo "script usage: $(basename \$0) [-o ParentOneVCF] [-t ParentTwoVCF] [-c ParentOneChainFile] [-h ParentTwoChainFile] [-r repeatsBEDfile]  [-S SNPsTableOut]"
    esac
done

# if script arguemnets are not provided, exit
if [ $# -eq 0 ]; then
    >&2 echo "No arguments provided"
    exit 1
fi

# SNPs_input <- "/oasis/tscc/scratch/mioverto/LOH_methods/Overton_etal/variants/reconciled/all_BYxRM_clean_SNPs.RData"
# train_ID <- "00"
# fractionTrainIDs <- 1
# fractionNonTrainIDs <- 0.5
# output_dir <- "/oasis/tscc/scratch/mioverto/LOH_methods/Overton_etal/variants/reconciled/"
# output_format <- "R"

# End Section #################################################################
###############################################################################

module load R


reconcile_script=${codeDir}/hybridCall/05b_reconcile_callsets.R
# reconcile_script=${codeDir}/hybridCall/05b_reconcile_callsets_build.R

# P1_VCF=${finalVarDir}/${P1}/all_${P1}.sif.vcf
# P2_VCF=${finalVarDir}/${P2}/all_${P2}.sif.vcf
# P1_chain=${refDir}/final/${P1}/${P1}x${typeRef}.chain
# P2_chain=${refDir}/final/${P2}/${P2}x${typeRef}.chain
# repeats_bed=${refDir}/final/${typeRef}/repeats/${typeRef}_repeats.bed
# # SNPs_output=${refDir}/final/${typeRef}/${P1}_${P2}_FinalSNPs.RData
# SNPs_output=${finalVarDir}/${P1}_${P2}_SNPs.csv

Rscript $reconcile_script -o $P1_VCF -t $P2_VCF -c $P1_chain -h $P2_chain -r $repeats_bed -S $SNPs_output