#!/bin/bash
#PBS -A mioverto
#PBS -l nodes=1 
#PBS -l walltime=24:00:00

while getopts 'A:C:O:' OPTION; do
    case "$OPTION" in
    A)
        acc="$OPTARG"
        echo "Accession list is $OPTARG"
        ;;
    C)
        cache="$OPTARG"
        echo "SRA prefetch cache $OPTARG"
        ;;
    O)
        outDir="$OPTARG"
        echo "Output reads to $OPTARG"
        ;;
    ?)
        echo "script usage: $(basename \$0) [-A accessions] [-C prefetchCache] [-O outputDirectory]"
    esac
done


# pCache=$fasterqCache
# outDir=$readsDir
# accID=ERR5324367
# proj=test
# export preDir=/oasis/tscc/scratch/mioverto/LOH_methods/${proj}/reads/pre
# export outDir=/oasis/tscc/scratch/mioverto/LOH_methods/${proj}/reads/raw

# cd $preDir
# $prefetch $accID -O $preDir/$accID
# $fasterq -S $accID -O $outDir

cd $cache

if [[ -f ${acc} ]]; then
    accBase=$(basename "${acc}")
    echo "downloading reads from accessions in ${accBase}"
    while IFS= read -r accNum || [[ -n "$accNum" ]]; do
        export accID=$accNum    
        if [[ -f ${outDir}/${accID}_1.fastq ]]; then
            echo "Read files already exist"
        else
            echo "downloading ${accID}"
            $prefetch -O ${cache} ${accID} 
            $fasterq -S ${accID} -O ${outDir}
            rm  -rd ${homeCache}/${accID}
            rm  -rd ${cache}/${accID}
        fi
    done < "$acc"

    else
    echo "downloading reads for accession $acc"
            $prefetch -O $cache $acc 
            $fasterq -S $acc -O $outDir
            rm  -rd ${homeCache}/${acc}
            rm  -rd ${cache}/${acc}
fi
