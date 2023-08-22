

proj=Dutta_etal_2021

# export rawDir=/oasis/tscc/scratch/mioverto/geneDrive/reads/raw
#${seqRun}
# readsTrimDir=${readsRawDir/raw/test}
# export trimDir=/oasis/tscc/scratch/mioverto/geneDrive/reads/trim
export trimDir=/oasis/tscc/scratch/mioverto/LOH_methods/${proj}/parents/reads/raw

export logDir=/oasis/tscc/scratch/mioverto/LOH_methods/log/trim
export qcDir=/oasis/tscc/scratch/mioverto/LOH_methods/${proj}/parents/reads/QC

module load fastqc
# trimDir=${pReadsDir}
# qcDir=${pReadsDir}
for rfile in ${trimDir}/ABA_R*; do
    echo $rfile
    fastqc -o ${qcDir} -extract $rfile 
done


```
for r1file in ${trimDir}/*_R1*; do
    #i=$(($i+1))
    export R1COMP=$r1file
    export R2COMP=${R1COMP/R1/R2}
    export tmp=$(basename "${R1COMP/_R1/}")
    export index=${tmp:0:3}
    echo submitting ${index} 
# done
    qsub \
        -V \
        -o ${logDir}/trim_${index}_${DATE}.out \
        -e ${logDir}/trim_${index}_${DATE}.err \
        -N trim_${index} \
        ${script}
done


fromDir=/oasis/tscc/scratch/mioverto/geneDrive/reads/QC

toDir=./
inc="F_A00_*"

rclone copy TSCC:$fromDir --include $inc $toDir

```