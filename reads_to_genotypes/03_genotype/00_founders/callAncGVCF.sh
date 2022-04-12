#!/bin/bash
#PBS -A mioverto
#PBS -l nodes=1 
#PBS -l walltime=40:00:00

PATH=/usr/lib/jvm/java-1.8.0-openjdk-1.8.0.252.b09-2.el7_8.x86_64/bin:$PATH
GATK=/home/mioverto/bin/gatk_4/gatk-package-4.1.8.0-local.jar

# Call variants with HaplotypeCaller. -ERC generates a gVCF file
java -jar $GATK HaplotypeCaller  \
    -R ${refSeq} \
    -I ${bamRG} \
    -O ${VCFout} \
    -ERC GVCF 