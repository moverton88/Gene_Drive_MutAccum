'''
We require a alternative reference sequence for alignment of reads to both 
parental genomes and avoid mapping bias. 
'''

module load bwa
module load samtools
module load bcftools
module load bedtools

# BY reference fasta
oldRefSeq=/home/mioverto/geneDrive/refseq/BY/S288C_R64_refseq.fna
# Query RM scaffolds
altSeq=/home/mioverto/geneDrive/refseq/RM/other_refs/RM11-1a_composite.fna
#altSeq=/home/mioverto/geneDrive/refseq/RM/UCI_scaffolds/RM11-1a_UCI_2019.fna

# Alignment output
bamFile=/oasis/tscc/scratch/mioverto/geneDrive/refseq/RM/RM_comp_refseq.bam
bamCrct=${bamFile/.bam/.crct.bam}

VCFout=/home/mioverto/geneDrive/POS_files/RMxBY_comp_ref_HC.vcf
VCFnrm=${VCFout/bcf/nrm}

altRefSeq=/home/mioverto/geneDrive/refseq/RM/RM_refseq_comp_2021.fna

#Switch to Java 1.8 for GATK v4.X to function correctly. The long dir name may have to be updated when java is updated on the cluster
PATH=/usr/lib/jvm/java-1.8.0-openjdk-1.8.0.252.b09-2.el7_8.x86_64/bin:$PATH
# Path to GATK application
GATK=/home/mioverto/bin/gatk_4/gatk-package-4.1.8.0-local.jar

#************************************************************************
# Align BY reference with RM scaffolds
# Convert sam to bam, sort by positions, and index
# Call variants with bcftools pipeline - with indels called
bwa mem $oldRefSeq $altSeq | samtools view -S -b -r RMxBY | samtools sort -o $bamFile 
samtools index $bamFile
# bcftools mpileup --open-prob 1 -e 1 --excl-flags DUP -f $oldRefSeq $BAMCRCT > ${VCFOUT/.vcf/.test.vcf}

java -jar $GATK AddOrReplaceReadGroups \
   -I $bamFile \
   -O $bamCrct \
   -R $oldRefSeq \
   -LB lib1 \
   -PL ILLUMINA \
   -PU assembly \
   -SM RMxBY

samtools index $bamCrct

bcftools mpileup -f $oldRefSeq $bamCrct | bcftools call -m -v -Ov -o $VCFout

bcftools norm -m +any -f $oldRefSeq $VCFout -Ov -o ${VCFnrm}

bgzip -c ${VCFnrm} > $VCFgz
bcftools index $VCFgz
bcftools consensus -H A -f $oldRefSeq $VCFgz > $altRefSeq
# bcftools consensus -H A -f $oldRefSeq -c $chainOut $VCFgz > $altRefSeq

```
# If bam file lacks required metadata, this should help to troubleshoot
java -jar $GATK ValidateSamFile \
     -I $SAMFILE \
     -R $REFIN
    
java -jar $GATK AddOrReplaceReadGroups \
   -I $BAMFILE \
   -O $BAMCRCT \
   -R $REFIN \
   -LB lib1 \
   -PL ILLUMINA \
   -PU assembly \
   -SM RMxBY

samtools index $BAMCRCT

```
#######################################################################

```
REFDIC=/home/mioverto/geneDrive/refseq/BY/S288C_R64_refseq.dict

# GATK requires its own indexing prep
java -jar $GATK CreateSequenceDictionary \
    -R $REFIN \
    -O $REFDIC


# Call variants between BY and RM
java -jar $GATK HaplotypeCaller  \
   --read-filter AllowAllReadsReadFilter \
   --disable-read-filter WellformedReadFilter \
   -R $REFIN \
   -I $BAMCRCT \
   -O $VCFHC



java -jar $GATK IndexFeatureFile \
   -I $VCFOUT

# Generate an RM reference sequence from BY reference and VCF
java -jar $GATK FastaAlternateReferenceMaker  \
   -R $REFIN \
   -O $REFOUT \
   -V $VCFOUT


```



```
REFOUT=/oasis/tscc/scratch/mioverto/data/refseq/RM_ref/RM_refseq_UCSD_2020_v3.fna
rclone copy TSCC:$REFOUT ./

rclone copy TSCC:/oasis/tscc/scratch/mioverto/data/refseq/RM_ref/bam/vcf/RMxBY_ref_bcf.vcf ./
```