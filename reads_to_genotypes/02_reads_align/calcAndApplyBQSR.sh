#!/bin/bash
#PBS -A mioverto
#PBS -l nodes=1 
#PBS -l walltime=5:00:00


PATH=/usr/lib/jvm/java-1.8.0-openjdk-1.8.0.252.b09-2.el7_8.x86_64/bin:$PATH
GATK=/home/mioverto/bin/gatk_4/gatk-package-4.1.8.0-local.jar


# If bam file lacks required metadata, this should help to troubleshoot
java -Xmx16G -jar $GATK ValidateSamFile \
     -I ${bamDeDup} \
     -R ${refSeq}
    
java -Xmx16G -jar $GATK AddOrReplaceReadGroups \
   -I ${bamDeDup} \
   -O ${bamCrct} \
   -R ${refSeq} \
   -LB lib1 \
   -PL ILLUMINA \
   -PU assembly \
   -SM ${index}


# samtools index $BAMCRCT

java -Xmx16G -jar $GATK BaseRecalibrator \
   -I ${bamCrct} \
   -R ${refSeq} \
   --known-sites ${vcfIn} \
   -O ${bamTable}

java -Xmx16G -jar $GATK ApplyBQSR \
   -R ${refSeq} \
   -I ${bamCrct} \
   --bqsr-recal-file ${bamTable} \
   -O ${bamOut}