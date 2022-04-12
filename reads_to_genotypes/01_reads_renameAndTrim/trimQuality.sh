

export rawDir=/oasis/tscc/scratch/mioverto/geneDrive/reads/raw
#${seqRun}
# readsTrimDir=${readsRawDir/raw/test}
export trimDir=/oasis/tscc/scratch/mioverto/geneDrive/reads/trim
export logDir=/oasis/tscc/scratch/mioverto/geneDrive/log/trim
export qcDir=/oasis/tscc/scratch/mioverto/geneDrive/reads/QC

module load fastqc

for r1file in ${trimDir}/N_A00*R1*; do
    echo $r1file
    fastqc -o ${qcDir} -extract $r1file 
done
    #i=$(($i+1))
    export R1COMP=$r1file
    export R2COMP=${R1COMP/R1/R2}
    export tmp=$(basename "${R1COMP/_R1/}")
    export index=${tmp:0:5}
    echo submitting ${index} 
# done
    qsub \
        -V \
        -o ${logDir}/trim_${index}_${DATE}.out \
        -e ${logDir}/trim_${index}_${DATE}.err \
        -N trim_${index} \
        ${script}
done

```
fromDir=/oasis/tscc/scratch/mioverto/geneDrive/reads/QC

toDir=./
inc="F_A00_*"

rclone copy TSCC:$fromDir --include $inc $toDir

```