

PATH=/usr/lib/jvm/java-1.8.0-openjdk-1.8.0.252.b09-2.el7_8.x86_64/bin:$PATH
GATK=/home/mioverto/bin/gatk_4/gatk-package-4.1.8.0-local.jar


alignRef=RM
# export refseq=/home/mioverto/geneDrive/refseq/RM/RM_refseq_UCSD_2020_v4.fna
# export refseq=/home/mioverto/geneDrive/refseq/BY/S288C_R64_refseq.fna
export VCFdir=/oasis/tscc/scratch/mioverto/geneDrive/dualRef/${alignRef}_aligned/variants/allVarVcfs
export vcfFile=${VCFdir}/A_${alignRef}a${alignRef}c_allVar.vcf
export indelOut=${vcfFile/allVar.vcf/indels.vcf}
export SNPOut=${vcfFile/allVar.vcf/SNPs.vcf}

java -jar -Xmx16g $GATK SelectVariants  \
    -R $refseq \
    --variant $vcfFile \
    --select-type-to-include INDEL \
    -O ${indelOut}
    
java -jar -Xmx16g $GATK SelectVariants  \
    -R $refseq \
    --variant $vcfFile \
    --select-type-to-include SNP \
    -O ${SNPOut}