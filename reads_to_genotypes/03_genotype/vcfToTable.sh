#!/bin/bash
#PBS -A mioverto
#PBS -l nodes=1 
#PBS -l walltime=40:00:00

PATH=/usr/lib/jvm/java-1.8.0-openjdk-1.8.0.252.b09-2.el7_8.x86_64/bin:$PATH
GATK=/home/mioverto/bin/gatk_4/gatk-package-4.1.8.0-local.jar


alignRef=RM
VCFname=A_${alignRef}a${alignRef}c_allVar.vcf
VCFdir=/oasis/tscc/scratch/mioverto/geneDrive/dualRef/${alignRef}_aligned/variants/allVarVcfs
VCFin=${VCFdir}/${VCFname}
tableDir=${VCFdir}/tables
tableOut=${tableDir}/${VCFname/.vcf/.tsv}


# Call variants with HaplotypeCaller. -ERC generates a gVCF file
java -jar $GATK VariantsToTable \
     -V $VCFin \
     -F CHROM -F POS -F REF -F ALT -F QUAL \
     -F BaseQRankSum -F MQRankSum -F MQ -F QD -F FS -F SOR \
     -GF AD -GF GT -GF GQ \
     -O $tableOut