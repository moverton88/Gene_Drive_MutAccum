##################################################################
##################################################################
# CHOOSE REFERENCE SEQUENCE ######################################
# RM = RM reference
# BY = BY reference
# BYm = BY masked reference
# Cas9 = Cas9 gRNA construct
export alignRef=RM
export callRef=${alignRef}

export alignTag="default"
export tag="allVar"
# _local

# export intervalFile=/home/mioverto/geneDrive/POS_files/${callRef}_POS.bed


if [[ ${callRef} == RM ]] ; then
        export refSeq=/home/mioverto/geneDrive/refseq/RM/RM_refseq_UCSD_2020_v4.fna
        echo "found RM ref"
    elif [[ ${callRef} == BY ]] ; then
        export refSeq=/home/mioverto/geneDrive/refseq/BY/S288C_R64_refseq.fna
        echo "found BY ref"
    else
    echo "alignment reference not found"
fi

projDir=/oasis/tscc/scratch/mioverto/LOH_methods/Pankajam_etal_2020
# projDir=/oasis/tscc/scratch/mioverto/geneDrive/dualRef
if [[ ${alignRef} == BY || ${alignRef} == RM ]] ; then
        export bamDir=${projDir}/${alignRef}_aligned/bam
        export gVCFdir=${projDir}/${alignRef}_aligned/variants/gVCFs
        export lineVCFdir=${projDir}/${alignRef}_aligned/variants/gVCFs/lineMulti
        export finalVCFdir=${projDir}/${alignRef}_aligned/variants/allVarVcfs
    elif [[ ${alignRef} == Cas9_N ]] ; then
        export bamDir=/oasis/tscc/scratch/mioverto/geneDrive/Cas9/bam
        export gVCFdir=/oasis/tscc/scratch/mioverto/geneDrive/Cas9/variants/gVCFs
        export refSeq=/home/mioverto/geneDrive/refseq/Drive/GFP_NrsR_cassette.fna
        export lineVCFdir=/oasis/tscc/scratch/mioverto/geneDrive/Cas9/variants/gVCFs/lineMulti
        export finalVCFdir=/oasis/tscc/scratch/mioverto/geneDrive/Cas9/variants/allVarVcfs
    elif [[ ${alignRef} == Cas9_H ]] ; then
        export bamDir=/oasis/tscc/scratch/mioverto/geneDrive/Cas9/bam
        export gVCFdir=/oasis/tscc/scratch/mioverto/geneDrive/Cas9/variants/gVCFs
        export refSeq=/home/mioverto/geneDrive/refseq/Drive/Cas9_GFP_NrsR.fna
        export lineVCFdir=/oasis/tscc/scratch/mioverto/geneDrive/Cas9/variants/gVCFs/lineMulti
        export finalVCFdir=/oasis/tscc/scratch/mioverto/geneDrive/Cas9/variants/allVarVcfs
    elif [[ ${alignRef} == Cas9_F ]] ; then
        export bamDir=/oasis/tscc/scratch/mioverto/geneDrive/Cas9/bam/
        export gVCFdir=/oasis/tscc/scratch/mioverto/geneDrive/Cas9/variants/gVCFs/
        export refSeq=/home/mioverto/geneDrive/refseq/Drive/gRNA_Cas9_GFP_NrsR_WTade2.fna
        export lineVCFdir=/oasis/tscc/scratch/mioverto/geneDrive/Cas9/variants/gVCFs/lineMulti
        export finalVCFdir=/oasis/tscc/scratch/mioverto/geneDrive/Cas9/variants/allVarVcfs
    else
        echo "alignment reference does not exist"
fi


export script=/home/mioverto/code/variants/02_combineAndCall_gatkGgVCFs.sh
export logDir=/oasis/tscc/scratch/mioverto/LOH_methods/log/
# export logDir=/oasis/tscc/scratch/mioverto/geneDrive/log/callGvcf
export DATE=$(date +'%m_%d_%Y')

# Submitting jobs in a loop for files that have not been created yet
# Wildcard must include only founder "00" ID, as the script bases 
# lineage groups on this
lineage=L

for gVCF in ${gVCFdir}/*.g.vcf; do
    export index=$(basename "${gVCF}" .g.vcf)
    # export lineage=${index:0:3}
    # export VCFout=${finalVCFdir}/${lineage}_${alignRef}a${callRef}c_${alignTag}_${tag}.vcf
    export VCFout=${finalVCFdir}/${lineage}_${alignRef}_${alignTag}_${tag}.vcf
    echo Submitting combine and call $(basename "${VCFout}")
# done
    qsub \
        -V \
        -N genotype_${lineage} \
        -o ${logDir}/multi-gVCF_${index}_${DATE}.out \
        -e ${logDir}/multi-gVCF_${index}_${DATE}.err \
        ${script}
done


```
export tag="_allVar"
export refseq=/home/mioverto/geneDrive/refseq/Cas9_gRNA_construct.fasta
export gVCFdir=/oasis/tscc/scratch/mioverto/geneDrive/Cas9/variants/gVCFs/
export lineVCFdir=/oasis/tscc/scratch/mioverto/geneDrive/Cas9/variants/gVCFs/
export finalVCFdir=/oasis/tscc/scratch/mioverto/geneDrive/Cas9/variants/finalVCFs
export VCFout=${finalVCFdir}/${lineage}_${alignRef}_${tag}.vcf




Missing anc:
*N_B00 - fastq misnamed N_B01 - fixed
N_F00 - seq failed
*N_G00 .idx - redo
N_H00 - seq failed

Poor anc:
N_C00
N_E00
H_F00
H_H00
F_A00

Poor data:
N_A01 - small fastq
N_A04 - nrml fastq
N_B09 - small fastq
N_C03
N_D07
N_E01
N_E11
N_F04
N_F07
N_F09
N_H02
N_H04
H_A12 - good fastq
H_B09 - good fastq
H_B10 - good fastq
H_B12 - good fastq
H_C06
H_C07
H_C09
H_D03
H_E01
H_E03
H_E11
H_F02
H_F04
H_F07
H_F12
H_G07
H_H03
F_A07 - good fastq
F_B01 - good fastq
F_B06 - good fastq
F_B07 - good fastq
F_B14 - good fastq
F_C06
F_D05
F_E03
F_E06
F_E09
F_F09 *

```