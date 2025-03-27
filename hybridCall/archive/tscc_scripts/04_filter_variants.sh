#!/bin/bash
#PBS -A mioverto
#PBS -l nodes=1 
#PBS -l walltime=5:00:00

while getopts 'v:l:f:o:R:' OPTION; do
   case "$OPTION" in
    v)
        variantFile="$OPTARG"
        echo "Variant file $OPTARG"
    ;;

   R)
      refSeq="${OPTARG}"
      echo "Operating system $OPTARG"
      ;;

   ?)
        echo "script usage: $(basename \$0) [-v variantFile] [-R referenceSequence]"
    esac
done

# lineage=${line}
# refSeq=$P1refSeq
refFile=$(basename "${refSeq}")
refName=${refFile%%_*}

# variantFile=${finalVarDir}/${refName}/${lineage}_${refName}${tag}.final.vcf
# variantFile=${variantsDir}/individuals/${refName}/L_${P1}${tag}.vcf
# variantFile=${finalVCF}
variantPrefix=${variantFile/.final.vcf/}
variantOut=${variantFile/.final.vcf/.filter.vcf}
# variantPrefix=${variantFile/.vcf/}
# variantOut=${variantPrefix}.filter.vcf

java -jar -Xmx16g $gatk SelectVariants  \
    -R ${refSeq} \
    --variant ${variantFile} \
    --select-type-to-include INDEL \
    -O ${variantPrefix}.indel.vcf

java -jar -Xmx16g $gatk SelectVariants  \
    -R ${refSeq} \
    --variant ${variantFile} \
    --select-type-to-include SNP \
    -O ${variantPrefix}.snps.vcf

# rm ${variant}

# SNP filter
java -jar $gatk VariantFiltration \
    -R $refSeq \
    -V ${variantPrefix}.snps.vcf \
    -O ${variantPrefix}.snps.m.vcf \
    --filter-name "QD2" --filter-expression "QD < 2.0" \
    --filter-name "FS60" --filter-expression "FS > 60.0" \
    --filter-name "MQ20" --filter-expression "MQ < 20.0" \
    --filter-name "QUAL30" --filter-expression "QUAL < 30.0" \
    --filter-name "SOR10" --filter-expression "SOR > 10.0" \
    --filter-name "DP6" --filter-expression "DP < 6"

java -jar $gatk SelectVariants \
    -V ${variantPrefix}.snps.m.vcf \
    -O ${variantPrefix}.snps.f.vcf \
    --exclude-filtered true

# Indel filter
java -jar $gatk VariantFiltration \
    -R $refSeq \
    -V ${variantPrefix}.indel.vcf \
    -O ${variantPrefix}.indel.m.vcf \
    --filter-name "QD2" --filter-expression "QD < 2.0" \
    --filter-name "FS200" --filter-expression "FS > 200.0" \
    --filter-name "QUAL30" --filter-expression "QUAL < 30.0" \
    --filter-name "DP6" --filter-expression "DP < 6"

java -jar $gatk SelectVariants \
    -V ${variantPrefix}.indel.m.vcf \
    -O ${variantPrefix}.indel.f.vcf \
    --exclude-filtered true 

java -jar -Xmx16g $gatk MergeVcfs  \
   -I ${variantPrefix}.snps.f.vcf \
   -I ${variantPrefix}.indel.f.vcf  \
   -O ${variantOut}

rm ${variantPrefix}.indel.vcf*
rm ${variantPrefix}.indel.m.vcf*
rm ${variantPrefix}.indel.f.vcf*

rm ${variantPrefix}.snps.vcf*
rm ${variantPrefix}.snps.m.vcf*
rm ${variantPrefix}.snps.f.vcf*
