#!/bin/bash

# Global variables ############################################################
export opSys=tscc
export baseDir=/oasis/tscc/scratch/mioverto/LOH_methods
export appDir=/home/mioverto/bin
export codeDir=/home/mioverto/code/hybridCall

## Overton et al 2022
export proj=Overton_etal
export projDir=${baseDir}/${proj}
export typeRef=S288C
export P1=BY
export P2=RM

## Sui et al 2020
export proj=Sui_etal
export projDir=${baseDir}/${proj}
export typeRef=S288C
export P1=W303
export P2=YJM789
export line=L

## Dutta et al 2021
export proj=Dutta_etal
export projDir=${baseDir}/${proj}
export typeRef=S288C

export P1=ABS
export P2=BKL
export line=H01

export P1=ABP
export P2=BFQ
export line=H02

export P1=BAP
export P2=BAN
export line=H03

export P1=BTI
export P2=ABA
export line=H04

export P1=ACD
export P2=AKQ
export line=H05

export P1=ACK
export P2=CMQ
export line=H06

export P1=ACG
export P2=BAK
export line=H07

export P1=CGD
export P2=AKE
export line=H08

export P1=BAM
export P2=CPG
export line=H09

export set_vars=/home/mioverto/code/hybridCall/00b_set_variables.sh
source $set_vars ${projDir} ${typeRef} ${P1} ${P2}

# Run 01_initialize_project.sh ################################################

export ini_proj=/home/mioverto/code/hybridCall/00a_initialize_project.sh

# chmod +x $ini_proj

if [[ ! -d ${baseDir}/${proj} ]]; then
    echo "creating project directories and files"
    export ini_proj=~/home/mioverto/code/hybridCall/01_initialize_project.sh
    $ini_proj \
        -s $opSys \
        -p $proj \
        -b $baseDir \
        -a $appDir \
        -c $codeDir \
        -t $typeRef \
        -O $P1 \
        -T $P2
    else
    echo "directory found"
fi


# End #########################################################################
###############################################################################
#                                                                             #
#                                                                             #
###############################################################################
# Run 00b_set_variables.sh ####################################################

export set_vars=/home/mioverto/code/hybridCall/00b_set_variables.sh
# chmod +x $set_vars

source $set_vars ${projDir} ${typeRef} ${P1} ${P2}
# source $set_vars ${projDir} S288C ${P2} ${P1}

# End #########################################################################
###############################################################################
#                                                                             #
#                                                                             #
###############################################################################
# Run 01_download_SRA_reads.sh ################################################

export download_SRA=/home/mioverto/code/hybridCall/01_download_SRA_reads.sh
# chmod +x $download_SRA

# Parent sequences
# pAcc=ERR1309072
pAcc=${metaDir}/Dutta_etal_2021_parent_accessions.txt
qsub -V -N ${proj}_download_sra_parent \
    -o  ${projDir}/logs/sra_$(date +'%Y_%m_%d').out \
    -e  ${projDir}/logs/sra_$(date +'%Y_%m_%d').err \
    -F "-A ${pAcc} -C ${fasterqCache} -O ${pReadsDir}" \
    $download_SRA

# Sample sequences
# acc=${metaDir}/Sui_etal_2020_miss_accessions.txt
qsub -V -N ${proj}_download_sra_clones \
    -o  ${projDir}/logs/sra_$(date +'%Y_%m_%d').out \
    -e  ${projDir}/logs/sra_$(date +'%Y_%m_%d').err \
    -F "-A ${acc} -C ${fasterqCache} -O ${readsDir}" \
    $download_SRA

# End #########################################################################
###############################################################################
#                                                                             #
#                                                                             #
###############################################################################
# Run rename_fastq.sh ##################################################

export rename=/home/mioverto/code/hybridCall/rename_fastq.sh
# chmod +x $rename

# export plan_file=${metaDir}/Sui_etal_2020_planfile_crct.csv
export plan_file=${metaDir}/Dutta_etal_2021_parent_rename.csv

qsub -V -N rename_files \
-o  ${projDir}/logs/rename_files_$(date +'%Y_%m_%d').out \
-e  ${projDir}/logs/rename_files_$(date +'%Y_%m_%d').err \
-F  "-r ${readsDir}/ -p ${plan_file}" \
${rename}

${rename} -r ${pReadsDir} -p ${plan_file}
# End #########################################################################
###############################################################################
#                                                                             #
#                                                                             #
###############################################################################
# Run trim_reads.sh #########################################################

# export trim_reads=/home/mioverto/LOH_methods/code/trim_reads.sh
# # chmod +x $trim_reads

# export plan_file=${metaDir}/Dutta_etal_2021_EPclones_rename.csv

# qsub -V -N rename_files \
# -o  ${projDir}/logs/rename_files_$(date +'%Y_%m_%d').out \
# -e  ${projDir}/logs/rename_files_$(date +'%Y_%m_%d').err \
# -F  "-r ${readsDir} -p ${plan_file}" \
# ${trim_reads}

# End #########################################################################
###############################################################################
#                                                                             #
#                                                                             #
###############################################################################
# Run construct_refSeq.sh #####################################################

export construct_refSeq=/home/mioverto/code/hybridCall/construct_refSeq.sh
# export construct_refSeq=/home/mioverto/LOH_methods/code/construct_refSeq_contigs.sh
# chmod +x $construct_refSeq

P_in=${P1}
export pReads=${pReadsDir}/${P_in}_R1.fastq
# export pReads=${pReadsDir}/${P_in}/${P_in}_contigs.fasta
export typeRefSeq=${typeRefDir}/S288C_R64_refseq.fasta

qsub -V -N ${P_in}_construct_refseq \
        -o  ${projDir}/logs/${P_in}_constRef_$(date +'%Y_%m_%d').out \
        -e  ${projDir}/logs/${P_in}_constRef_$(date +'%Y_%m_%d').err \
        -F  "-r ${pReads} -T ${typeRefSeq} -s tscc" \
        $construct_refSeq

# $construct_refSeq -r ${pReads} -T ${typeRefSeq} -s tscc

# End #########################################################################
###############################################################################
#                                                                             #
#                                                                             #
###############################################################################
# Run 02_align_and_call.sh ####################################################

export align_call=/home/mioverto/code/hybridCall/02_align_and_call.sh
# chmod +x $align_call

line_P1_P2=${projDir}/metadata/hybrids_and_parents.csv

while IFS="," read -r line P1 P2 || [[ -n "$line" ]]; do

    export $line
    export $P1
    export $P2

    source $set_vars ${projDir} ${typeRef} ${P1} ${P2}

    for readsFile in ${readsDir}/${line}*R1*fastq; do
        # export readsFile=${readsDir}/H02_09_R1.fastq
        export readsName=$(basename "${readsFile%%_R1*}")
        export sampleName=${readsName}_${P1}${tag}
        if [ -f ${alignDir}/${P1}/${readsName}*.dedup.bam ]; then
            echo "${sampleName} alignment already exists, skipping"

        else
            echo "Generating ${sampleName} alignment"
            
            qsub -V -N ${sampleName}_align \
                -o  ${projDir}/logs/${sampleName}_align_$(date +'%Y_%m_%d').out \
                -e  ${projDir}/logs/${sampleName}_align_$(date +'%Y_%m_%d').err \
                -F  "-p ${P1} -r ${readsFile} -v ${indvVarDir}/${P1} -t gVCF -R ${P1refSeq} -s tscc" \
                ${align_call}
        fi
    done
 
    for readsFile in ${readsDir}/${line}*R1*fastq; do
        # export readsFile=${readsDir}/H02_09_R1.fastq
        export readsName=$(basename "${readsFile%%_R1*}")
        export sampleName=${readsName}_${P2}${tag}
        if [ -f ${alignDir}/${P2}/${readsName}*.dedup.bam ]; then
            echo "${sampleName} alignment already exists, skipping"

        else
            echo "Generating ${sampleName} alignment"
            
            qsub -V -N ${sampleName}_align \
                -o  ${projDir}/logs/${sampleName}_align_$(date +'%Y_%m_%d').out \
                -e  ${projDir}/logs/${sampleName}_align_$(date +'%Y_%m_%d').err \
                -F  "-p ${P2} -r ${readsFile} -v ${indvVarDir}/${P2} -t gVCF -R ${P2refSeq} -s tscc" \
                ${align_call}
        fi
    done
done < ${line_P1_P2}

# End #########################################################################
###############################################################################
#                                                                             #
#                                                                             #
###############################################################################
# Run 03_combineAndCall_gVCFs.sh ####################################################

export combine_call=/home/mioverto/code/hybridCall/03_combineAndCall_gVCFs.sh
# chmod +x $combine_call

line_P1_P2=${projDir}/metadata/hybrids_and_parents.csv

while IFS="," read -r line P1 P2 || [[ -n "$line" ]]; do
    echo ${line} ${P1} ${P2}
    source $set_vars ${projDir} ${typeRef} ${P1} ${P2}
    export sampleName=${line}_${P1}

    if [ -f ${finalVarDir}/${P1}/${sampleName}*.final.vcf ]; then
        echo "${sampleName} alignment already exists, skipping"

    else
        echo "Generating ${sampleName} alignment"

        qsub -V -N ${sampleName}_combine_call \
            -o  ${projDir}/logs/${sampleName}_combineCall_$(date +'%Y_%m_%d').out \
            -e  ${projDir}/logs/${sampleName}_combineCall_$(date +'%Y_%m_%d').err \
            -F  "-v ${indvVarDir}/${P1} -l ${line} -f "none" -p ${P1} -R ${P1refSeq}" \
            ${combine_call}
    fi

    export sampleName=${line}_${P2}

    if [ -f ${finalVarDir}/${P2}/${sampleName}*.final.vcf ]; then
        echo "${sampleName} alignment already exists, skipping"

    else
        echo "Generating ${sampleName} alignment"
  
        qsub -V -N ${sampleName}_combine_call \
            -o  ${projDir}/logs/${sampleName}_combineCall_$(date +'%Y_%m_%d').out \
            -e  ${projDir}/logs/${sampleName}_combineCall_$(date +'%Y_%m_%d').err \
            -F  "-v ${indvVarDir}/${P2} -l ${line} -f "none" -p ${P2} -R ${P2refSeq}" \
            ${combine_call}

    fi

done < ${line_P1_P2}

# End #########################################################################
###############################################################################
#                                                                             #
#                                                                             #
###############################################################################
# Run 04_filter_variants.sh ####################################################

export filter_variants=/home/mioverto/code/hybridCall/04_filter_variants.sh
# chmod +x $filter_variants

line_P1_P2=${projDir}/metadata/hybrids_and_parents.csv

while IFS="," read -r line P1 P2 || [[ -n "$line" ]]; do

    source $set_vars ${projDir} ${typeRef} ${P1} ${P2}

    sampleName=${line}_${P1}
    VCFin=${finalVarDir}/${P1}/${sampleName}.final.vcf
    qsub -V -N ${sampleName}_filter_variants \
                -o  ${projDir}/logs/${sampleName}_filterVariants_$(date +'%Y_%m_%d').out \
                -e  ${projDir}/logs/${sampleName}_filterVariants_$(date +'%Y_%m_%d').err \
                -F  "-v ${VCFin} -R ${P1refSeq}" \
                ${filter_variants}

    sampleName=${line}_${P2}
    VCFin=${finalVarDir}/${P2}/${sampleName}.final.vcf
    qsub -V -N ${sampleName}_filter_variants \
                -o  ${projDir}/logs/${sampleName}_filterVariants_$(date +'%Y_%m_%d').out \
                -e  ${projDir}/logs/${sampleName}_filterVariants_$(date +'%Y_%m_%d').err \
                -F  "-v ${VCFin} -R ${P2refSeq}" \
                ${filter_variants}

done < ${line_P1_P2}

# End #########################################################################
###############################################################################
#                                                                             #
#                                                                             #
###############################################################################
# Run allele_depths_atBEDintervals.sh ####################################################

export allele_DP=/home/mioverto/code/hybridCall/allele_depths_atBEDintervals.sh
# chmod +x $allele_DP

export set_vars=/home/mioverto/code/hybridCall/00b_set_variables.sh
source $set_vars ${projDir} ${typeRef} ${P1} ${P2}
# source $set_vars ${projDir} ${typeRef} ${P2} ${P1}

export IDs_list=${metaDir}/AD_IDs.tsv
export bed_file=${metaDir}/hom_sites.bed
# export line=N_A

 while IFS= read -r ID || [[ -n "$ID" ]]; do
        # export accID=$accNum
        bam_file=${alignDir}/${P1}/${ID}_${P1}.dedup.bam
        tsv_out=${variantsDir}/final/${P1}/AD_tables/${ID}_${P1}_homAD.tsv
        # echo ${bam_file}
        # echo ${tsv_out}
        if [[ ! -d ${variantsDir}/final/${P1}/AD_tables ]]; then
            echo "Creating table directory"
            mkdir ${variantsDir}/final/${P1}/AD_tables
        fi
        if [[ -f ${tsv_out} ]]; then
            echo "Allele depth table for $(basename ${tsv_out%_${P1}*}) already exist"
        else
            qsub -V -N $(basename ${tsv_out%.tsv}) \
                -o  ${projDir}/logs/${tsv_out/.tsv/}_$(date +'%Y_%m_%d').out \
                -e  ${projDir}/logs/${tsv_out/.tsv/}_$(date +'%Y_%m_%d').err \
                -F  "-p ${P1} -a ${bam_file} -b ${bed_file} -o ${tsv_out} -R ${P1refSeq}" \
                ${allele_DP}
        fi
    done < "${IDs_list}"


# ID=N_C06
# bam_file=${alignDir}/${P1}/${ID}_${P1}.dedup.bam
# tsv_out=${variantsDir}/final/${P1}/${ID}_${P1}_homAD.tsv
# ${allele_DP} -p ${P1} -a ${bam_file} -b ${bed_file} -o ${tsv_out} -R ${P1refSeq}

# End #########################################################################
###############################################################################
#                                                                             #
#                                                                             #
###############################################################################
# Run 05_reconcile_callsets.sh ####################################################

export reconcile_calls=/home/mioverto/code/hybridCall/05_reconcile_callsets.sh
# chmod +x $reconcile_calls

export set_vars=/home/mioverto/code/hybridCall/00b_set_variables.sh

line_P1_P2=${projDir}/metadata/hybrids_and_parents.csv

while IFS="," read -r line P1 P2 || [[ -n "$line" ]]; do

    source $set_vars ${projDir} ${typeRef} ${P1} ${P2}
    P1_VCF=${finalVarDir}/${P1}/${line}_${P1}.filter.vcf
    P2_VCF=${finalVarDir}/${P2}/${line}_${P2}.filter.vcf
    P1_chain=${refDir}/final/${P1}/${P1}x${typeRef}.chain
    P2_chain=${refDir}/final/${P2}/${P2}x${typeRef}.chain
    repeats_bed=${refDir}/final/${typeRef}/repeats/${typeRef}_repeats.bed
    SNPs_output=${variantsDir}/reconciled/${line}_${P1}x${P2}_reconciled_SNPs.RData

qsub -V -N ${line}_${P1}x${P2}_reconcile_calls \
                -o  ${projDir}/logs/reconcile/reconcile_${P1}x${P2}_$(date +'%Y_%m_%d').out \
                -e  ${projDir}/logs/reconcile/reconcile_${P1}x${P2}_$(date +'%Y_%m_%d').err \
                -F  "-o $P1_VCF -t $P2_VCF -c $P1_chain -h $P2_chain -r $repeats_bed -S $SNPs_output" \
                ${reconcile_calls}

done < "${line_P1_P2}"

# $reconcile_calls -o $P1_VCF -t $P2_VCF -c $P1_chain -h $P2_chain -r $repeats_bed -S $SNPs_output

# End #########################################################################
###############################################################################
#                                                                             #
#                                                                             #
###############################################################################
# Run remove_prob_samples.sh ####################################################

export remove_prob_samples=/home/mioverto/code/hybridCall/remove_prob_samples.sh
# chmod +x $remove_prob_samples

SNPs_input=${variantsDir}/reconciled/${samples}_${P1}x${P2}_reconciled_SNPs.RData
SNPs_output=${variantsDir}/reconciled/${samples}_${P1}x${P2}_clean_SNPs.RData

module load R
Rscript $remove_prob_samples -s $SNPs_input -S $SNPs_output

# End #########################################################################
###############################################################################
#                                                                             #
#                                                                             #
###############################################################################
# Run 06a_train_het_model.R ####################################################

# set_vars=/home/mioverto/code/hybridCall/00b_set_variables.sh
source $set_vars ${projDir} ${typeRef} ${P1} ${P2}

export train_hets=/home/mioverto/code/hybridCall/06a_train_het_model.R
# chmod +x $train_hets

module load R

line_P1_P2=${projDir}/metadata/hybrids_and_parents.csv

while IFS="," read -r line P1 P2 || [[ -n "$line" ]]; do
    echo $line
    # SNPs_input=${variantsDir}/reconciled/${line}_${P1}x${P2}_clean_SNPs.RData
    SNPs_input=${variantsDir}/reconciled/${line}_${P1}x${P2}_reconciled_SNPs.RData
    IDs=$line
    fract_train=0.75
    fract_nontrain=0
    input_params="4,1000,100"
    model_output=${variantsDir}/reconciled/

    # qsub -V -N train_hets \
    #                 -o  ${projDir}/logs/reconcile/reconcile_${P1}x${P2}_$(date +'%Y_%m_%d').out \
    #                 -e  ${projDir}/logs/reconcile/reconcile_${P1}x${P2}_$(date +'%Y_%m_%d').err \
    #                 -F  "-s $SNPs_input -I $IDs -f $fract_train -F $fract_nontrain -p $input_params -o $model_output -O R" \
    #                 ${train_hets}

    Rscript $train_hets -s $SNPs_input -I $IDs -f $fract_train -F $fract_nontrain -p $input_params -o $model_output -O R
done < ${line_P1_P2}

# End #########################################################################
###############################################################################
#                                                                             #
#                                                                             #
###############################################################################
# Run 06b_filter_het_calls.R ####################################################

filter_hets=/home/mioverto/code/hybridCall/06b_filter_het_calls.R
export filter_het_calls=/home/mioverto/code/hybridCall/06_filter_het_calls.sh
# chmod +x $filter_het_calls
set_vars=/home/mioverto/code/hybridCall/00b_set_variables.sh
model_input=/oasis/tscc/scratch/mioverto/LOH_methods/Overton_etal/variants/reconciled/all_BYxRM_hetTrain.RData
# model_input=/oasis/tscc/scratch/mioverto/LOH_methods/Sui_etal/variants/reconciled/L_W303xYJM789_hetTrain.RData

line_P1_P2=${metaDir}/hybrids_and_parents.csv
while IFS="," read -r line P1 P2 || [[ -n "$line" ]]; do
    echo $line
    source $set_vars ${projDir} ${typeRef} ${P1} ${P2}
    samples=${line}
    SNPs_input=${variantsDir}/reconciled/${samples}_${P1}x${P2}_reconciled_SNPs.RData
    # model_input=${variantsDir}/reconciled/all_BYxRM_hetTrain.RData
    model_type="v"
    FE_rate=0.05
    SNPs_output=${SNPs_input}

    qsub -V -N filter_hets \
                    -o  ${projDir}/logs/reconcile/filter_hets_${P1}x${P2}_$(date +'%Y_%m_%d').out \
                    -e  ${projDir}/logs/reconcile/filter_hets_${P1}x${P2}_$(date +'%Y_%m_%d').err \
                    -F  "-s $SNPs_input -m $model_input -t $model_type -r $FE_rate -S $SNPs_output" \
                    ${filter_het_calls}

done < ${line_P1_P2}

module load R
Rscript $filter_hets -s $SNPs_input -m $model_input -t $model_type -r $FE_rate -S $SNPs_output

# End #########################################################################
###############################################################################
#                                                                             #
#                                                                             #
###############################################################################
# Run 06a_train_hom_model.R ####################################################

export train_homs=/home/mioverto/code/hybridCall/06a_train_hom_model.R
# chmod +x $train_homs

module load R

samples="all"
hom_AD_dir=${variantsDir}/hom_AD/
# SNPs_input=${variantsDir}/reconciled/${P1}x${P2}_clean_SNPs.csv
IDs="00"
fract_train=0.9
fract_nontrain=0.5
model_output=${variantsDir}/reconciled/

# qsub -V -N train_hets \
#                 -o  ${projDir}/logs/reconcile/reconcile_${P1}x${P2}_$(date +'%Y_%m_%d').out \
#                 -e  ${projDir}/logs/reconcile/reconcile_${P1}x${P2}_$(date +'%Y_%m_%d').err \
#                 -F  "-s $SNPs_input -I $IDs -f $fract_train -F $fract_nontrain -o $model_output -O R" \
#                 ${train_homs}

Rscript $train_homs -h $hom_AD_dir -c $IDs -fn $fract_train -O $fract_nontrain


# End #########################################################################
###############################################################################
#                                                                             #
#                                                                             #
###############################################################################
# Run 06b_filter_hom_calls.R ####################################################

export filter_homs=/home/mioverto/code/hybridCall/06b_filter_hom_calls.R
export filter_hom_calls=/home/mioverto/code/hybridCall/06b_filter_hom_calls.sh
# chmod +x $filter_homs
# chmod +x $filter_hom_calls


# line="H01"

# SNPs_input=${variantsDir}/reconciled/${P1}x${P2}_clean_SNPs.csv
model_input=/oasis/tscc/scratch/mioverto/LOH_methods/Overton_etal/variants/reconciled/hom_model.RData
FE_rate=0.05

line_P1_P2=${metaDir}/hybrids_and_parents.csv
while IFS="," read -r line P1 P2 || [[ -n "$line" ]]; do
    echo $line
    source $set_vars ${projDir} ${typeRef} ${P1} ${P2}
    samples=${line}
    SNPs_input=${variantsDir}/reconciled/${samples}_${P1}x${P2}_reconciled_SNPs.RData
    SNPs_output=${SNPs_input}

    qsub -V -N filter_homs \
                    -o  ${projDir}/logs/reconcile/homFilter_${line}_${P1}x${P2}_$(date +'%Y_%m_%d').out \
                    -e  ${projDir}/logs/reconcile/homFilter_${line}_${P1}x${P2}_$(date +'%Y_%m_%d').err \
                    -F  "-s $SNPs_input -I $IDs -f $fract_train -F $fract_nontrain -o $model_output -O R" \
                    ${filter_homs}

done < ${line_P1_P2}

module load R

SNPs_input=${variantsDir}/reconciled/${line}_${P1}x${P2}_reconciled_SNPs.RData
SNPs_output=${SNPs_input}

Rscript $filter_homs -s $SNPs_input -m $model_input -r $FE_rate -S $SNPs_output

