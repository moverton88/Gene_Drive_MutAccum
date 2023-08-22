
export refSeq=/home/mioverto/geneDrive/refseq/BY/S288C_R64_refseq.fna

export bamDir=/oasis/tscc/scratch/mioverto/geneDrive/dualRef/BY_aligned/bam/DeDup
export bamFile=${bamDir}/anc_BYc_local.sort.rg.bam

export VCFdir=/oasis/tscc/scratch/mioverto/geneDrive/dualRef/BY_aligned/variants/allVarVcfs/
export VCFout=${VCFdir}/anc_pool_BYc_local.vcf

export script=/home/mioverto/code/variants/call_AltHC.sh
export logDir=/oasis/tscc/scratch/mioverto/geneDrive/log/callGvcf
export DATE=$(date +'%m_%d_%Y')

qsub \
    -V \
    -N alt_call \
    -o ${logDir}/multi-gVCF_${index}_${DATE}.out \
    -e ${logDir}/multi-gVCF_${index}_${DATE}.err \
    ${script}