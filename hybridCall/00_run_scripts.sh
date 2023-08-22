#!/bin/bash

# Global variables ############################################################
export opSys=tscc
export baseDir=/oasis/tscc/scratch/mioverto/LOH_methods
export appDir=/home/mioverto/bin
export codeDir=/home/mioverto/code


## Overton et al 2022
export proj=Overton_etal
export projDir=${baseDir}/${proj}
export typeRef=S288C
export P1=BY
export P2=RM


## Sui et al 2020
export proj=Sui_etal_2020
export projDir=${baseDir}/${proj}
export typeRef=S288C
export P1=W303
export P2=YJM789

## test
export proj=test
export projDir=${baseDir}/${proj}
export typeRef=S288C

# export P1=ABS
# export P2=BKL
# export line=H01

# export P1=ABP
# export P2=BFQ
# export line=H02

export P1=BAP
export P2=BAN
export line=H03

# export P1=BTI
# export P2=ABA
# export line=H04

# export P1=ACD
# export P2=AKQ
# export line=H05

# export P1=ACK
# export P2=CMQ
# export line=H06

# export P1=ACG
# export P2=BAK
# export line=H07

# export P1=CGD
# export P2=AKE
# export line=H08

# export P1=BAM
# export P2=CPG
# export line=H09


```
# Run 01_initialize_project.sh ################################################

export ini_proj=/home/mioverto/LOH_methods/code/00a_initialize_project.sh

# chmod +x $ini_proj

if [[ ! -d ${baseDir}/${proj} ]]; then
    echo "creating project directories and files"
    export ini_proj=~/LOH_methods/code/01_initialize_project.sh
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


```

###############################################################################
# Run 02_load_variables.sh ####################################################

export set_vars=/home/mioverto/LOH_methods/code/00b_set_variables.sh
# chmod +x $set_vars

source $set_vars ${projDir} S288C ${P1} ${P2}
# source $set_vars ${projDir} S288C ${P2} ${P1}

# End #########################################################################
###############################################################################
#                                                                             #
#                                                                             #
###############################################################################
# Run 01_download_SRA_reads.sh ################################################

export download_SRA=/home/mioverto/LOH_methods/code/01_download_SRA_reads.sh
# chmod +x $download_SRA

qsub -V -N ${proj}_download_sra \
    -o  ${projDir}/logs/sra_$(date +'%Y_%m_%d').out \
    -e  ${projDir}/logs/sra_$(date +'%Y_%m_%d').err \
    -F "-L ${accFile} -C ${fasterqCache} -O ${readsDir}" \
    $download_SRA

# End #########################################################################
###############################################################################
#                                                                             #
#                                                                             #
###############################################################################
# Run 03_construct_refSeq.sh ##################################################

export construct_refSeq=/home/mioverto/LOH_methods/code/03_construct_refSeq.sh
# export construct_refSeq=/home/mioverto/LOH_methods/code/03_construct_refSeq_contigs.sh
# chmod +x $construct_refSeq

export set_vars=/home/mioverto/LOH_methods/code/00b_set_variables.sh
source $set_vars ${projDir} ${typeRef} ${P1} ${P2}
# source $set_vars ${projDir} ${typeRef} ${P2} ${P1}

for parent in ${pReadsDir}/${P1}*_R1.fastq; do
    export tmp=$(basename "${parent}" .fastq)
    export pName=${tmp%%_*}
    if [ -d "${refSeqDir}/${pName}" ]; then
        echo "${pName} parental reference already exists, skipping"
    else
        export P1=${pName}
        P2=X
        source $set_vars ${projDir} ${typeRef} ${P1} ${P2}
        echo "Generating ${pName} parental reference"
        
        qsub -V -N ${P1}_construct_refseq \
        -o  ${projDir}/logs/${P1}_constRef_$(date +'%Y_%m_%d').out \
        -e  ${projDir}/logs/${P1}_constRef_$(date +'%Y_%m_%d').err \
        $construct_refSeq
    fi
done

# End #########################################################################
###############################################################################
#                                                                             #
#                                                                             #
###############################################################################
# Run 04_align_markDups.sh ####################################################

export align_markDups=/home/mioverto/LOH_methods/code/04_align_markDups.sh
# chmod +x $align_markDups

# export line=N_A00
# export tag=_mum
# export readsFile=${readsDir}/F_A00_R1.trim.fastq
# export readsName=$(basename "${readsFile}" _R1.trim.fastq)
# export sampleName=${readsName}_${P1}${tag}
# export P1=S288C


# export lineFile=${metaDir}/hybrids_and_parents.csv
# export lineFile=${metaDir}/Sui_clone_IDs.txt
# export lineFile=${metaDir}/expIDs_lineageIDs.N.csv

# grep "N_" ${lineFile} > ${lineFile/.csv/.N.csv}

export set_vars=/home/mioverto/LOH_methods/code/00b_set_variables.sh
source $set_vars ${projDir} ${typeRef} ${P1} ${P2}
# P1refSeq=${P1refSeq/.fna/${_tag}.fna}
# P1=S288C
# P1refSeq=${typeRefSeq}
# while IFS="," read -r line P1 P2 || [[ -n "$line" ]]; do
# while IFS="," read -r line || [[ -n "$line" ]]; do
while IFS="," read -r old line || [[ -n "$line" ]]; do
    # source $set_vars ${projDir} ${typeRef} ${P1} ${P2}
    echo $line $P1 $P2

    for readsFile in ${readsDir}/${line}*R1.trim.fastq; do
        # export readsFile=${readsDir}/N_A00_R1P.trim.fastq
        export readsName=$(basename "${readsFile%%_R1*}")
        export sampleName=${readsName}_${P1}${tag}
        if [ -f "${alignDir}/${P1}/${sampleName}.dedup.bam" ]; then
            echo "${sampleName} alignment already exists, skipping"

        else
            echo "Generating ${sampleName} alignment"
            
            qsub -V -N ${sampleName}_align \
            -o  ${projDir}/logs/${sampleName}_align_$(date +'%Y_%m_%d').out \
            -e  ${projDir}/logs/${sampleName}_align_$(date +'%Y_%m_%d').err \
            -F  "-r ${readsFile} -a ${alignDir}/${P1} -R ${P1refSeq} -s set_vars" \
            ${align_markDups}
        fi
    done

done < "$lineFile"

source $set_vars ${projDir} ${typeRef} ${P2} ${P1}

# while IFS="," read -r line P1 P2 || [[ -n "$line" ]]; do
while IFS="," read -r line || [[ -n "$line" ]]; do
    for readsFile in ${readsDir}/${line}*R1.trim.fastq; do
        export readsName=$(basename "${readsFile%%_R1*}")
        export sampleName=${readsName}_${P1}${tag}
        if [ -f "${alignDir}/${P1}/${sampleName}.dedup.bam" ]; then
            echo "${sampleName} alignment already exists, skipping"

        else
            echo "Generating ${sampleName} alignment"
            
            qsub -V -N ${sampleName}_align \
            -o  ${projDir}/logs/${sampleName}_align_$(date +'%Y_%m_%d').out \
            -e  ${projDir}/logs/${sampleName}_align_$(date +'%Y_%m_%d').err \
            -F  "-r ${readsFile} -a ${alignDir}/${P1} -R ${P1refSeq} -s set_vars" \
            ${align_markDups}
        fi
    done
done < "$lineFile"

# End #########################################################################
###############################################################################
#                                                                             #
#                                                                             #
###############################################################################
# Run 05_call_gVcf_gatkHC.sh ##################################################

export call_gVcf=/home/mioverto/LOH_methods/code/05_call_gVcf_gatkHC.sh
# chmod +x $call_gVcf

# P1=ABS
export set_vars=/home/mioverto/LOH_methods/code/00b_set_variables.sh
source $set_vars ${projDir} ${typeRef} ${P1} ${P2}

export lineFile=${metaDir}/hybrids_and_parents_3.csv
# export line=H
# export tag="_MQ30"
# export tag="_mum"
# P1refSeq=${P1refSeq/.fna/_mum.fna}

while IFS="," read -r line P1 P2 || [[ -n "$line" ]]; do
    source $set_vars ${projDir} ${typeRef} ${P1} ${P2}
    echo $line $P1 $P2

    for alignFile in ${alignDir}/${P1}/${line}*${P1}${tag}.dedup.bam; do
        # alignFile=${alignDir}/${P1}/H03_03_${P1}${tag}.dedup.bam
        export sampleName=$(basename "${alignFile}" .dedup.bam)
        export alignName=${sampleName/_${P1}${tag}/}
        echo "looking for ${sampleName} variants"
        if [ -f "${variantsDir}/individuals/${P1}/${sampleName}.g.vcf" ]; then
            echo "${sampleName} gVCF already exists, skipping"

        else
            echo "Generating ${sampleName} gVCF"
        
            qsub -V -N ${sampleName}_call_gVCF \
            -o  ${projDir}/logs/${sampleName}_callgVCF_$(date +'%Y_%m_%d').out \
            -e  ${projDir}/logs/${sampleName}_callgVCF_$(date +'%Y_%m_%d').err \
            -F  "-a ${alignFile} -v ${variantsDir}/individuals/${P1} -R ${P1refSeq} -t g" \
            ${call_gVcf}

        fi

    done

    source $set_vars ${projDir} ${typeRef} ${P2} ${P1}

    for sample in ${alignDir}/${P1}/${line}*.dedup.bam; do
        export sampleName=$(basename "${sample}" .dedup.bam)

        if [ -f "${variantsDir}/individuals/${P1}/${sampleName}${tag}.g.vcf" ]; then
            echo "${sampleName} gVCF already exists, skipping"

        else
            echo "Generating ${sampleName} gVCF"
        
            qsub -V -N ${sampleName}_call_gVCF \
            -o  ${projDir}/logs/${sampleName}_${P1}_callgVCF_$(date +'%Y_%m_%d').out \
            -e  ${projDir}/logs/${sampleName}_${P1}_callgVCF_$(date +'%Y_%m_%d').err \
            -F  "-a ${sample} -R ${P1refSeq} -s set_vars" \
            ${call_gVcf}

        fi
        
    done
done < "$lineFile"

# End #########################################################################
###############################################################################
#                                                                             #
#                                                                             #
###############################################################################
# Run 05b_filterArtifacts.sh ##################################################

export filter_artifacts=/home/mioverto/LOH_methods/code/05b_filterArtifacts.sh
# chmod +x $filter_artifacts
# export P1=ABA
# export set_vars=/home/mioverto/LOH_methods/code/00b_set_variables.sh
# source $set_vars ${projDir} ${typeRef} ${P1} ${P2}
# export sampleName=L000
# export line=L

# export lineFile=${metaDir}/hybrids_and_parents.csv
export lineFile=${metaDir}/Sui_clone_IDs.txt
# P1=ABS
export set_vars=/home/mioverto/LOH_methods/code/00b_set_variables.sh
source $set_vars ${projDir} ${typeRef} ${P1} ${P2}

export lineFile=${metaDir}/hybrids_and_parents_3.csv
export set_vars=/home/mioverto/LOH_methods/code/00b_set_variables.sh

while IFS="," read -r line P1 P2 || [[ -n "$line" ]]; do
    source $set_vars ${projDir} ${typeRef} ${P1} ${P2}
    echo $line $P1 $P2

    for sample in ${alignDir}/${P1}/${line}*.g.vcf; do
        # sample=${variantsDir}/${P1}/${line}_01_${P1}.g.vcf
        # sample=${variantsDir}/individuals/${P1}/${sampleName}_${P1}.g.vcf
        alignFile=${alignDir}/${P1}/${sampleName}_${P1}.dedup.bam
        export sampleName=$(basename "${sample}" .g.vcf)
        # echo "looking for ${sampleName} variants"
        if [ -f "${variantsDir}/individuals/${P1}/${sampleName}.faa.g.vcf" ]; then
            echo "${sampleName} FAA gVCF already exists, skipping"

        else
            echo "Generating ${sampleName} FAA gVCF"
        
            qsub -V -N ${sampleName}_FAA \
            -o  ${projDir}/logs/${sampleName}_FAA_$(date +'%Y_%m_%d').out \
            -e  ${projDir}/logs/${sampleName}_FAA_$(date +'%Y_%m_%d').err \
            -F  "-v ${sample} -a ${alignFile} -R ${P1refSeq} -s set_vars" \
            ${filter_artifacts}

        fi

    done

    source $set_vars ${projDir} ${typeRef} ${P2} ${P1}
    echo $line $P1 $P2

    for sample in ${alignDir}/${P1}/${line}*.g.vcf; do
        # sample=${variantsDir}/${P1}/${line}_01_${P1}.g.vcf
        # sample=${variantsDir}/individuals/${P1}/${sampleName}_${P1}.g.vcf
        export sampleName=$(basename "${sample}" .g.vcf)
        # echo "looking for ${sampleName} variants"
        if [ -f "${variantsDir}/individuals/${P1}/${sampleName}.faa.g.vcf" ]; then
            echo "${sampleName} FAA gVCF already exists, skipping"

        else
            echo "Generating ${sampleName} FAA gVCF"
        
            qsub -V -N ${sampleName}_FAA \
            -o  ${projDir}/logs/${sampleName}_FAA_$(date +'%Y_%m_%d').out \
            -e  ${projDir}/logs/${sampleName}_FAA_$(date +'%Y_%m_%d').err \
            -F  "-a ${sample} -R ${P1refSeq} -s set_vars" \
            ${filter_artifacts}

        fi

    done
done < "$lineFile"

# End #########################################################################
###############################################################################
#                                                                             #
#                                                                             #
###############################################################################
# Run 06_combineAndCall_gVCFs.sh ##############################################

export cmbnAndCall=/home/mioverto/LOH_methods/code/06_combineAndCall_gVCFs.sh
# chmod +x $cmbnAndCall

export set_vars=/home/mioverto/LOH_methods/code/00b_set_variables.sh
source $set_vars ${projDir} ${typeRef} ${P1} ${P2}

# export founders=${metaDir}/founderList.txt
# tag=_MQ30
# line=all
# ${cmbnAndCall} -v ${variantsDir}/individuals/${P1} -l ${line} -f 'none' -R ${P1refSeq}

sampleName=${line}_${P1}${tag}
qsub -V -N ${line}_${P1}${tag}_cmbnAndCall \
-o  ${projDir}/logs/${sampleName}_cmbnAndCall_$(date +'%Y_%m_%d').out \
-e  ${projDir}/logs/${sampleName}_cmbnAndCall_$(date +'%Y_%m_%d').err \
-F  "-v ${variantsDir}/individuals/${P1} -l ${line} -f $founders -R ${P1refSeq}" \
${cmbnAndCall}

source $set_vars ${projDir} ${typeRef} ${P2} ${P1}

# ${cmbnAndCall} -v ${variantsDir}/individuals/${P1} -l ${line} -f 'none' -R ${P1refSeq}

qsub -V -N ${line}_${P1}_cmbnAndCall \
-o  ${projDir}/logs/${sampleName}_cmbnAndCall_$(date +'%Y_%m_%d').out \
-e  ${projDir}/logs/${sampleName}_cmbnAndCall_$(date +'%Y_%m_%d').err \
-F  "-v ${variantsDir}/individuals/${P1} -l ${line} -f $founders -R ${P1refSeq}" \
${cmbnAndCall}
