#!/bin/bash
#PBS -A mioverto
#PBS -l nodes=1 
#PBS -l walltime=5:00:00


module load bowtie2
module load samtools
module load bcftools

PATH=/usr/lib/jvm/java-1.8.0-openjdk-1.8.0.252.b09-2.el7_8.x86_64/bin:$PATH
GATK=/home/mioverto/bin/gatk_4/gatk-package-4.1.8.0-local.jar

# Start with alignment. Assumes renamed, trimmed fastq files and un-indexed reference fasta file

# Reference must be indexed by bowtie2 before read alignment
if [ ! -f ${REFPREFIX}.1.bt2 ]; then
    echo "Reference sequence needs to be indexed (bowtie2-build)"
    bowtie2-build ${REFSEQ} ${REFPREFIX}
fi

# Align paired and unpaired reads reads with bowtie2. 
# ${INDEX} used for naming alignments. Output is sorted bam alignment
# -I, -X flag give the minimum, maximum valid gap between paired reads
# samtools view -h [include header] -b [bam output] -u [uncompressed .bam] | sort -m [memory allocated] -o [output file] -T [temp files]
bowtie2 --rg-id ${index} --rg SM:${index} --very-sensitive-local -I 100 -X 1000 \
-x ${REFPREFIX} -1 ${R1PFILE} -2 ${R2PFILE} -U "${R1UFILE},${R2UFILE}" \
| samtools view -h -b -u | samtools sort -m 10000000 -o ${bamRaw} -O bam -T ${bamDir}/temp_dir/${index}

# other bowtie2 options: -N 1 -L 20 --mp 3,2 --rdg 3,2 --rfg 3,2 --local
# Remove duplicate reads
java -jar $GATK MarkDuplicates \
	--INPUT ${bamRaw} \
	--OUTPUT ${bamDeDup} \
	--METRICS_FILE ${metrics}/${index}.txt \
	--ASSUME_SORTED TRUE \
    --TMP_DIR ${bamDir}/temp_dir/${index}

# index bam alignment
samtools index ${bamDeDup}

# rm -r ${bamDir}/temp_dir/${index}


```
java -Xmx16G -jar $GATK ValidateSamFile \
     -I ${bamRaw} \
     -R ${REFSEQ}

```