# Run script under interactive mode

export bamDir=/oasis/tscc/scratch/mioverto/geneDrive/dualRef/BY_aligned/bam/DeDup
export bamIn=${bamDir}/anc_BYc_local.bam
export bamOut=${bamDir}/anc_BYc_local.sort.bam
export bamRG=${bamDir}/anc_BYc_local.sort.rg.bam

anc_list=${bamDir}/anc_bam.list
dir ${bamDir}/*00*.bam >> ${anc_list}
# less ${anc_list}

module load samtools

samtools merge ${bamDir}/anc_BYc_local.bam -b ${anc_list}
samtools sort ${bamIn} -o ${bamOut}
samtools index -b ${bamOut} 

java -jar $GATK AddOrReplaceReadGroups \
    -I ${bamOut} \
    -O ${bamRG} \
    -SORT_ORDER coordinate \
    -RGID anc \
    -RGLB clones \
    -RGPL illumina \
    -RGPU lane.2 \
    -RGSM all \
    -CREATE_INDEX True