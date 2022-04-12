
export alignRef=BY
export alignRef=RM

export tag=local
# export alignTag=""
# _local


if [ ${alignRef} == BYm ] ; then
    export bamDir=/oasis/tscc/scratch/mioverto/geneDrive/ambiRef/bam/
    export gVCFdir=/oasis/tscc/scratch/mioverto/geneDrive/ambiRef/variants/gVCFs/${callRef}_call
    elif [[ ${alignRef} == BY || ${alignRef} == RM ]] ; then
        export bamDir=/oasis/tscc/scratch/mioverto/geneDrive/dualRef/${alignRef}_aligned/bam
        export gVCFdir=/oasis/tscc/scratch/mioverto/geneDrive/dualRef/${alignRef}_aligned/variants/gVCFs/${callRef}_call
    elif [[ ${alignRef} == Cas9 ]] ; then
        export bamDir=/oasis/tscc/scratch/mioverto/geneDrive/${alignRef}/bam
        export gVCFdir=/oasis/tscc/scratch/mioverto/geneDrive/${alignRef}/variants/gVCFs/
        export refSeq=/home/mioverto/geneDrive/refseq/Cas9_gRNA_construct.fasta
    else
        echo "alignment reference does not exist"
fi



if [ ${alignRef} != Cas9 ] ; then
    echo "looking for call ref ..."
    if [[ ${callRef} == RM ]] ; then
        export refSeq=/home/mioverto/geneDrive/refseq/RM/RM_refseq_UCSD_2020_v4.fna
        echo "found RM ref"
    elif [[ ${callRef} == BY ]] ; then
        export refSeq=/home/mioverto/geneDrive/refseq/BY/S288C_R64_refseq.fna
        echo "found BY ref"
    fi
    else
    callRef=${alignRef}
fi


export bamDir=/oasis/tscc/scratch/mioverto/geneDrive/dualRef/${alignRef}_aligned/bam/DeDup
export bamRG=${bamDir}/anc_${alignRef}c_${tag}.sort.rg.bam

export VCFdir=/oasis/tscc/scratch/mioverto/geneDrive/dualRef/${alignRef}_aligned/variants/allVarVcfs/
export VCFout=${VCFdir}/anc_pool_${alignRef}c_${tag}.vcf

export script=/home/mioverto/code/variants/callAncHC.sh
export logDir=/oasis/tscc/scratch/mioverto/geneDrive/log/callGvcf
export DATE=$(date +'%m_%d_%Y')

qsub \
    -V \
    -N anc_call \
    -o ${logDir}/multi-gVCF_${index}_${DATE}.out \
    -e ${logDir}/multi-gVCF_${index}_${DATE}.err \
    ${script}