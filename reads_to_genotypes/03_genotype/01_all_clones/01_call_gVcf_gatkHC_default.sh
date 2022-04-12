#!/bin/bash
#PBS -A mioverto
#PBS -l nodes=1 
#PBS -l walltime=7:00:00


PATH=/usr/lib/jvm/java-1.8.0-openjdk-1.8.0.252.b09-2.el7_8.x86_64/bin:$PATH
GATK=/home/mioverto/bin/gatk_4/gatk-package-4.1.8.0-local.jar

# Call variants with HaplotypeCaller. -ERC generates a gVCF file
if [ -z ${vcfIn} ]; then
    java  -Xmx16G -jar $GATK HaplotypeCaller \
        -R $refSeq \
        -I ${bamDeDup} \
        -O ${gVCFout} \
        -ERC GVCF 
fi

if [ ! -z ${vcfIn} ]; then
    java -Xmx16G -jar $GATK HaplotypeCaller \
        -R $refSeq \
        -I ${bamDeDup} \
        -alleles ${vcfIn} \
        -O ${gVCFout} \
        -ERC GVCF 
fi


``` 
bamOut=${gVCFdir}/realign_bam/${index}_${tag}${alignTag}_realign.bam
vcfIn=/home/mioverto/geneDrive/POS_files/RMvcf/RMxBY_ref_noMit.vcf
gVCFout=${gVCFdir}/${index}_local_allele.g.vcf

# Call variants with HaplotypeCaller. -ERC generates a gVCF file
java -jar $GATK HaplotypeCaller  \
    -R $refSeq \
    -I ${bamDeDup} \
    -bamout ${bamOut} \
    -O ${gVCFout} \
    --active-probability-threshold 0.001 \
    --min-dangling-branch-length 3 \
    --max-reads-per-alignment-start 100 \
    --max-assembly-region-size 500 \
    -ERC GVCF 



```

module load matlab/2021a

matlab ${simScript} -I $input -O $output

