#!/bin/bash

#SBATCH --partition=hotel
#SBATCH --QOS=hotel
#SBATCH --nodes=1
# ... Other SLURM options ...

# Global variables ############################################################
export opSys=slurm
export baseDir=/tscc/lustre/ddn/scratch/$USER/LOH_methods
export appDir=~/bin
export codeDir=~/code/hybridCall
export yyyymmdd=$(date +%F)

## Overton et al 2022
export proj=Overton_etal
export projDir=${baseDir}/${proj}
export typeRef=S288C
export P1=BY
export P2=RM
export line=all
export plan_fn=expIDs_lineageIDs.csv

## Sui et al 2020
export proj=Sui_etal
export projDir=${baseDir}/${proj}
export typeRef=S288C
export P1=W303
export P2=YJM789
export line=L
export accFile=Sui_etal_2020_accessions_crct.txt
export plan_fn=Sui_etal_2020_planfile_crct.csv

## Dutta et al 2021
export proj=Dutta_etal
export projDir=${baseDir}/${proj}
export typeRef=S288C
export accFile=Dutta_etal_2021_EPclonescrct_accessions_crct.txt
export pAccFile=Dutta_etal_2021_parent_accessions.txt
export plan_fn=Dutta_etal_2021_EPclones_rename.csv

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

## Pankajam al 2020
export proj=Pankajam_etal
export projDir=${baseDir}/${proj}
export typeRef=S288C
export P1=S288C
export P2=RM
export line=L
export accFile=Pankajam_etal_2020_accessions_SY.txt
export plan_fn=Pankajam_rename.csv

export P1=S288C
export P2=YJM789
export line=SY

## Loeillet et al 2020
export proj=Loeillet_etal
export projDir=${baseDir}/${proj}
export typeRef=S288C
export accFile=Loeillet_etal_2020_accessions_crct.txt
export plan_fn=Loeillet_rename.csv
export P1=BY
export P2=SK1
export line="all"

export set_vars=${HOME}/code/hybridCall/00b_set_variables_slurm.sh
source $set_vars ${projDir} ${typeRef} ${P1} ${P2}

# Run 01_initialize_project.sh ################################################

export ini_proj=${HOME}/code/hybridCall/00a_initialize_project.sh

# chmod +x $ini_proj

if [[ ! -d ${baseDir}/${proj} ]]; then
    echo "creating project directories and files"
    # export ini_proj=~/code/hybridCall/00a_initialize_project.sh
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
# Run 00b_set_variables_slurm.sh ##############################################

export set_vars=${HOME}/code/hybridCall/00b_set_variables_slurm.sh
# chmod +x $set_vars

source $set_vars ${projDir} ${typeRef} ${P1} ${P2}
# source $set_vars ${projDir} S288C ${P2} ${P1}

# End #########################################################################
###############################################################################
#                                                                             #
#                                                                             #
###############################################################################
# Run 01_download_SRA_reads.sh ################################################

export download_SRA=${HOME}/code/hybridCall/01_download_SRA_reads.sb
# chmod +x $download_SRA

# Parent sequences
# export acc=ERR475452
export acc=${metaDir}/${pAccFile}
export cache=${pFasterqCache}
export outDir=${pReadsDir}/${P2}
sbatch --export=ALL $download_SRA

# Sample sequences
export acc=${metaDir}/${accFile}
export cache=${fasterqCache}
export outDir=${readsDir}
sbatch --export=ALL $download_SRA

# End #########################################################################
###############################################################################
#                                                                             #
#                                                                             #
###############################################################################
# Run rename_fastq.sb ##################################################

export rename=${HOME}/code/hybridCall/rename_fastq.sb
# chmod +x $rename

export plan_file=${metaDir}/${plan_fn}

export rename_dir=${pReadsDir}
# export rename_dir=${readsDir}/raw
sbatch --export=ALL $rename


# ${rename} -r ${pReadsDir} -p ${plan_file}

# End #########################################################################
###############################################################################
#                                                                             #
#                                                                             #
###############################################################################
# Run QC_reads.sh #########################################################

export QC_reads=${HOME}/LOH_methods/code/QC_reads.sh
# chmod +x $QC_reads


# End #########################################################################
###############################################################################
#                                                                             #
#                                                                             #
###############################################################################
# Run trim_reads.sh #########################################################

export trim_reads=${HOME}/code/hybridCall/trim_reads.sh
# chmod +x $trim_reads

# reads1=${readsDir}/raw/F_A00_R1.fastq.gz  
export trimOpt=PE
export outputDir=${readsDir}
export crop=1000
export lead=10
export trail=10
export head=10
export adapters=/tscc/nfs/home/mioverto/bin/Trimmomatic-0.39/adapters/NexteraPE-PE.fa
for reads in ${readsDir}/raw/N*R1.fastq.gz; do
    # export reads=${readsDir}/raw/F_E00*R1.fastq.gz 
    export reads1=$reads
    echo $reads1
    export readsName=$(basename "${reads1%%_R1*}")
    if [ -f ${readsDir}/${readsName}*R1.fastq* ]; then
        echo "${readsName} already trimmed, skipping"

    else
        echo "Trimming reads from ${readsName}"
        
        sbatch --export=ALL $trim_reads
        # $trim_reads
    fi
done


# End #########################################################################
###############################################################################
#                                                                             #
#                                                                             #
###############################################################################
# Run construct_refSeq.sh #####################################################

export construct_refSeq=${HOME}/code/hybridCall/construct_refSeq.sb
# export construct_refSeq=${HOME}/code/hybridCall/construct_refSeq_contigs.sb
# chmod +x $construct_refSeq

P_in=${P2}
# P_in=${P2}
export pReads1=${pReadsDir}/${P_in}_R1.fastq
# export pReads=${pReadsDir}/${P_in}_contigs.fasta
export typeRefSeq=${typeRefDir}/S288C_refseq.fasta
# export typeRefSeq=${P1refSeq}
sbatch --export=ALL $construct_refSeq

# $construct_refSeq -r ${pReads} -T ${typeRefSeq} -s tscc

# End #########################################################################
###############################################################################
#                                                                             #
#                                                                             #
###############################################################################
# Run 02_align_and_call.sb ####################################################

export align_call=${HOME}/code/hybridCall/02_align_and_call.sb
# chmod +x $align_call

export P_in=${P1}
# export P_in=$typeRef
# export reads1=${readsDir}/wt_01_R1.fastq
export gVCFdir=${indvVarDir}/${P_in}
export VCFtype=gVCF
export refSeq=${refSeqDir}/${P_in}/${P_in}_refseq.fasta

# export line_P1_P2=${projDir}/metadata/hybrids_and_parents.csv

while IFS="," read -r line P1 P2 || [[ -n "$line" ]]; do

    source $set_vars ${projDir} ${typeRef} ${P1} ${P2}
    # line=wt
        for reads in ${readsDir}/N_*R1.*fastq*; do
            # export reads=${readsDir}/wt_00_R1.fastq
            export reads1=$reads
            echo $reads1
            export readsName=$(basename "${reads1%%_R1*}")
            export sampleName=${readsName}_${P_in}
            if [ -f ${alignDir}/${P_in}/${readsName}*.dedup.bam ]; then
                echo "${sampleName} alignment already exists, skipping"

            else
                echo "Generating ${sampleName} alignment"
                
                sbatch --export=ALL $align_call
                # $align_call
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
            
            sbatch --export=ALL $align_call --output=${projDir}/logs/alignCall_${proj}_${yyyymmdd}.out
            # $align_call
        fi
    done
done < ${line_P1_P2}

# End #########################################################################
###############################################################################
#                                                                             #
#                                                                             #
###############################################################################
# Run 03_combineAndCall_gVCFs.sh ##############################################

export combine_call=${HOME}/code/hybridCall/03_combineAndCall_gVCFs.sb
# chmod +x $combine_call

line_P1_P2=${projDir}/metadata/hybrids_and_parents.csv

export P_in=${P1}
# export line=all
export variantDir=${indvVarDir}/${P_in}
export founders="00"
export refSeq=${refSeqDir}/${P_in}/${P_in}_refseq.fasta

while IFS="," read -r line P1 P2 || [[ -n "$line" ]]; do
    echo ${line} ${P1} ${P2}
    source $set_vars ${projDir} ${typeRef} ${P1} ${P2}
    export sampleName=${line}_${P_in}

    if [ -f ${finalVarDir}/${P_in}/${line}_${P_in}.final.vcf ]; then
        echo "${sampleName} alignment already exists, skipping"

    else
        echo "Generating ${sampleName} alignment"

        sbatch --export=ALL $combine_call
            # $align_call
    fi

    export refSeq=${P2refSeq} 
    export variantDir=${indvVarDir}/${P2}
    export sampleName=${line}_${P2}

    if [ -f ${finalVarDir}/${P2}/${sampleName}*.final.vcf ]; then
        echo "${sampleName} alignment already exists, skipping"

    else
        echo "Generating ${sampleName} alignment"
  
        sbatch --export=ALL $combine_call

    fi

done < ${line_P1_P2}

# End #########################################################################
###############################################################################
#                                                                             #
#                                                                             #
###############################################################################
# Run 04_filter_variants.sb ####################################################

export filter_variants=${HOME}/code/hybridCall/04_filter_variants.sb
# chmod +x $filter_variants

line_P1_P2=${projDir}/metadata/hybrids_and_parents.csv


while IFS="," read -r line P1 P2 || [[ -n "$line" ]]; do

    source $set_vars ${projDir} ${typeRef} ${P1} ${P2}

    export refSeq=${P1refSeq}
    export variantFile=${finalVarDir}/${P1}/${line}_${P1}.final.vcf
    sbatch --export=ALL $filter_variants

    export refSeq=${P2refSeq}
    export variantFile=${finalVarDir}/${P2}/${line}_${P2}.final.vcf
    sbatch --export=ALL $filter_variants

done < ${line_P1_P2}

export P_in=$typeRef
export refSeq=${refSeqDir}/${P_in}/${P_in}_refseq.fasta
for file in ${finalVarDir}/${P_in}/*.final.vcf; do
    export line=${file%%_*}
    export variantFile=${file}
    if [ -f ${variantFile} ] && [ ! -f ${variantFile/.final./.filter.} ]; then
        echo "Filtering ${variantFile}"
        sbatch --export=ALL $filter_variants

    elif [ -f ${variantFile} ] && [ -f ${variantFile/.final./.filter.} ]; then
        echo "Filtered variant file already exists"
    else
        echo "variant file does not exist"
    fi
done


# End #########################################################################
###############################################################################
#                                                                             #
#                                                                             #
###############################################################################
# Run allele_depths_atBEDintervals.sh #########################################

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
# Run 05_00_VCFtoRData.sb ####################################################

# module spider r
# module load cpu/0.17.3 gcc/10.2.0-2ml3m2l r

export VCFtoRData=${HOME}/code/hybridCall/05_00_VCFtoRData.sb
# chmod +x $VCFtoRData
# line=all
line_P1_P2=${projDir}/metadata/hybrids_and_parents.csv

while IFS="," read -r line P1 P2 || [[ -n "$line" ]]; do

    source $set_vars ${projDir} ${typeRef} ${P1} ${P2}
    VCF_in=${finalVarDir}/${P1}/${line}_${P1}.filter.vcf
    RData_out=${VCF_in/filter.vcf/RData}
    if [ -f ${VCF_in} ] && [ ! -f ${RData_out} ]; then
        echo "reconcile single reference for file ${VCF_in} with line ${line}"
        sbatch --export=ALL $VCFtoRData

    elif [ -f ${VCF_in} ] && [ -f ${RData_out} ]; then
        echo "reconciled callset already exists for ${line}_${P1}"
        echo $line
    else
        echo "callset file does not exist"
    fi

    VCF_in=${finalVarDir}/${P2}/${line}_${P2}.filter.vcf
    RData_out=${VCF_in/filter.vcf/RData}

    if [ -f ${VCF_in} ] && [ ! -f ${RData_out} ]; then
        echo "reconcile single reference for file ${VCF_in} with line ${line}"
        sbatch --export=ALL $VCFtoRData

    elif [ -f ${VCF_in} ] && [ -f ${RData_out} ]; then
        echo "reconciled callset already exists for ${line}_${P2}"
        echo $line
    else
        echo "callset file does not exist"
    fi
 
done < "${line_P1_P2}"


# $reconcile_calls -o $P1_VCF -c $P1_chain -r $repeats_bed -S $SNPs_output

export P_in=$P2
# export P1_VCF=${finalVarDir}/${P_in}/${line}_${P_in}.filter.vcf
export P1_chain=${refDir}/final/${P_in}/${P_in}_chrom_lengths.tsv
# export P1_chain=${refDir}/final/${P_in}/${P_in}x${typeRef}.chain
export repeats_bed=${refDir}/final/${typeRef}/repeats/${typeRef}_repeats.bed
for file in ${finalVarDir}/${P_in}/*.filter.vcf; do
    export P1_VCF=$file
    P1_fn=$(basename ${file})
    export line=${P1_fn%%_*}
    export SNPs_output=${variantsDir}/final/${P_in}/${line}_${P_in}_singleRef_SNPs.RData
    if [ -f ${P1_VCF} ] && [ ! -f ${SNPs_output} ]; then
        echo "reconcile single reference for file ${P1_VCF} with line ${line}"
        sbatch --export=ALL $reconcile_calls

    elif [ -f ${P1_VCF} ] && [ -f ${SNPs_output} ]; then
        echo "reconciled callset already exists"
        echo $line
    else
        echo "callset file does not exist"
    fi
done

# End #########################################################################
###############################################################################
#                                                                             #
#                                                                             #
###############################################################################
# Run 05_reconcile_callsets.sb ####################################################

# module spider r
# module load cpu/0.17.3 gcc/10.2.0-2ml3m2l r

export reconcile_calls=${HOME}/code/hybridCall/05_reconcile_callsets.sb
# export reconcile_calls=${HOME}/code/hybridCall/05_reconcile_singleRef.sb
# chmod +x $reconcile_calls
# line=all
line_P1_P2=${projDir}/metadata/hybrids_and_parents.csv

while IFS="," read -r line P1 P2 || [[ -n "$line" ]]; do

    source $set_vars ${projDir} ${typeRef} ${P1} ${P2}
    export P1_VCF=${finalVarDir}/${P1}/${line}_${P1}.filter.vcf
    export P2_VCF=${finalVarDir}/${P2}/${line}_${P2}.filter.vcf
    export P1_chain=${refDir}/final/${P1}/${P1}x${typeRef}.chain
    export P2_chain=${refDir}/final/${P2}/${P2}x${typeRef}.chain
    export repeats_bed=${refDir}/final/${typeRef}/repeats/${typeRef}_repeats.bed
    export SNPs_output=${variantsDir}/final/Dual/${line}_${P1}x${P2}_reconciled_SNPs.RData

    sbatch --export=ALL $reconcile_calls
 
done < "${line_P1_P2}"


# $reconcile_calls -o $P1_VCF -c $P1_chain -r $repeats_bed -S $SNPs_output

export P_in=$P2
# export P1_VCF=${finalVarDir}/${P_in}/${line}_${P_in}.filter.vcf
export P1_chain=${refDir}/final/${P_in}/${P_in}_chrom_lengths.tsv
# export P1_chain=${refDir}/final/${P_in}/${P_in}x${typeRef}.chain
export repeats_bed=${refDir}/final/${typeRef}/repeats/${typeRef}_repeats.bed
for file in ${finalVarDir}/${P_in}/*.filter.vcf; do
    export P1_VCF=$file
    P1_fn=$(basename ${file})
    export line=${P1_fn%%_*}
    export SNPs_output=${variantsDir}/final/${P_in}/${line}_${P_in}_singleRef_SNPs.RData
    if [ -f ${P1_VCF} ] && [ ! -f ${SNPs_output} ]; then
        echo "reconcile single reference for file ${P1_VCF} with line ${line}"
        sbatch --export=ALL $reconcile_calls

    elif [ -f ${P1_VCF} ] && [ -f ${SNPs_output} ]; then
        echo "reconciled callset already exists"
        echo $line
    else
        echo "callset file does not exist"
    fi
done

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

# module load cpu/0.17.3 gcc/10.2.0-2ml3m2l r
Rscript $remove_prob_samples -s $SNPs_input -S $SNPs_output

# End #########################################################################
###############################################################################
#                                                                             #
#                                                                             #
###############################################################################
# Run 06a_train_het_model.R ####################################################

# set_vars=${HOME}/code/hybridCall/00b_set_variables_slurm.sh
source $set_vars ${projDir} ${typeRef} ${P1} ${P2}

export train_hets=/home/mioverto/code/hybridCall/06a_train_het_model.R
# chmod +x $train_hets

# module load cpu/0.17.3 gcc/10.2.0-2ml3m2l r

line_P1_P2=${projDir}/metadata/hybrids_and_parents.csv

while IFS="," read -r line P1 P2 || [[ -n "$line" ]]; do
    echo $line
    # SNPs_input=${variantsDir}/reconciled/${line}_${P1}x${P2}_clean_SNPs.RData
    SNPs_input=${variantsDir}/reconciled/${line}_${P1}x${P2}_reconciled_SNPs.RData
    IDs=$line
    fract_train=0.9 
    fract_nontrain=0.67
    input_params="4,1000,100"
    model_output=${variantsDir}/reconciled/

    # sbatch --export=ALL $XXXX

    Rscript $train_hets -s $SNPs_input -I $IDs -f $fract_train -F $fract_nontrain -p $input_params -o $model_output -O R
done < ${line_P1_P2}

# End #########################################################################
###############################################################################
#                                                                             #
#                                                                             #
###############################################################################
# Run 06b_filter_het_calls.R ####################################################

# export filter_hets=${HOME}/code/hybridCall/06b_filter_het_calls.R
export filter_het_calls=${HOME}/code/hybridCall/06b_filter_het_calls.sb
# chmod +x $filter_het_calls
set_vars=/home/mioverto/code/hybridCall/00b_set_variables.sh
model_input=/tscc/lustre/ddn/scratch/mioverto/LOH_methods/Overton_etal/variants/final/Dual/all_BYxRM_hetTrain.RData
# model_input=/tscc/lustre/ddn/scratch/mioverto/LOH_methods/Sui_etal/variants/final/Dual/L_W303xYJM789_hetTrain.RData
model_type="v"
FE_rate=0.05
    
line_P1_P2=${metaDir}/hybrids_and_parents.csv
while IFS="," read -r line P1 P2 || [[ -n "$line" ]]; do
    echo $line
    source $set_vars ${projDir} ${typeRef} ${P1} ${P2}
    SNPs_input=${variantsDir}/final/Dual/${line}_${P1}x${P2}_SNPs.RData
    # model_input=${variantsDir}/final/Dual/all_BYxRM_hetTrain.RData
    SNPs_output=${SNPs_input/_SNPs/_SNPs_hetFilter}

    # sbatch --export=ALL ${filter_het_calls}

done < ${line_P1_P2}

# module load cpu/0.17.3 gcc/10.2.0-2ml3m2l r
Rscript $filter_hets -s $SNPs_input -m $model_input -t $model_type -r $FE_rate -S $SNPs_output

export P_in=$typeRef
export model_input=/tscc/lustre/ddn/scratch/mioverto/LOH_methods/Overton_etal/variants/final/Dual/all_BYxRM_hetTrain.RData
export model_type="v"
export FE_rate=0.05
for file in ${finalVarDir}/${P_in}/*singlRef_SNPs.RData; do
    export line=${file%%_*}
    export SNPs_input=${finalVarDir}/${P_in}/${line}_${P_in}_singlRef_SNPs.RData
    # export SNPs_input=${file}
    export SNPs_output=${finalVarDir}/${P_in}/${line}_${P_in}_singlRef_SNPs_hf.RData
    # export SNPs_output=${file}
    if [ -f ${SNPs_input} ]; then
        echo "filter het calls in ${variantFile}"
        # sbatch --export=ALL $filter_het_calls
    else
        echo "callset file does not exist"
    fi
done


# End #########################################################################
###############################################################################
#                                                                             #
#                                                                             #
###############################################################################
# Run 06a_train_hom_model.R ####################################################

export train_homs=/home/mioverto/code/hybridCall/06a_train_hom_model.R
# chmod +x $train_homs

# module load cpu/0.17.3 gcc/10.2.0-2ml3m2l r

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
export filter_hom_calls=/home/mioverto/code/hybridCall/06_filter_hom_calls.sh
# chmod +x $filter_homs
# chmod +x $filter_hom_calls

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
                    -F  "-s $SNPs_input -m $model_input -r $FE_rate -S $SNPs_output" \
                    ${filter_hom_calls}

done < ${line_P1_P2}

# module load cpu/0.17.3 gcc/10.2.0-2ml3m2l r
# line="H01"
SNPs_input=${variantsDir}/reconciled/${line}_${P1}x${P2}_reconciled_SNPs.RData
SNPs_output=${SNPs_input}

Rscript $filter_homs -s $SNPs_input -m $model_input -r $FE_rate -S $SNPs_output

