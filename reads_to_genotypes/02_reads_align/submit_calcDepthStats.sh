

ref=RM
tag=sensitive

export bamDir=/oasis/tscc/scratch/mioverto/geneDrive/dualRef/${ref}_aligned/bam
export metricDir=${bamDir}/depth_metrics

if [[ ${ref} == RM ]] ; then
    export refSeq=/home/mioverto/geneDrive/refseq/RM/RM_refseq_UCSD_2020_v4.fna
    echo "found RM ref"
elif [[ ${ref} == BY ]] ; then
    export refSeq=/home/mioverto/geneDrive/refseq/BY/S288C_R64_refseq.fna
    echo "found BY ref"
fi

export script=/home/mioverto/code/align/calcDepthStats.sh

for bam in ${bamDir}/DeDup/*_local.dm.bam; do
    # export R1FILE=/oasis/tscc/scratch/mioverto/data/MAseq1/reads/trim/half-L100_1_R1P.trimmed.fastq
    export bamFile=${bam}
    export sampleName=$(basename "$bamFile" .dm.bam)
    export metricFile=${metricDir}/${sampleName}_wgs.txt
    echo "Submitting depth ${sampleName}"
# done
    qsub \
        -V \
        -N depth_${sampleName} \
        -o ${logDir}/depth_${sampleName}_${DATE}.out \
        -e ${logDir}/depth_${sampleName}_${DATE}.err \
        ${script}
done
