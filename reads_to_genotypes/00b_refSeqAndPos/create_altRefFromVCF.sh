'''
We require a RM reference sequence for alignment of reads to both parental genomes and avoid mapping bias. 

First, BWA aligns the scaffolds to the reference sequence and outputs a .sam file
We then sort and index the .sam file and output a .bam file
Bcftools mpileup to bcftools call creates the vcf file with all variant sites
between the reference sequence and alternative sequence
After normalizing the VCF format and indexing, bcftools consensus
creates and alternative reference sequence and bedtools maskfasta
replaces variant alleles with N characters
*  One thing to note is that I have only tested the masking function
   using a VCF that does not include indels. I do not know how this 
   function works on indels.
'''

module load bwa
module load samtools
module load bcftools
module load bedtools

# BY reference fasta
refIn=/home/mioverto/geneDrive/refseq/BY/S288C_R64_refseq.fna

# RM query sequence scaffolds
queryIn=/home/mioverto/geneDrive/refseq/RM/UCI_scaffolds/RM11-1a_UCI_2019.fna

# RM reference sequence output
newSeqOut=/home/mioverto/geneDrive/refseq/RM/RM_refseq_UCSD_2020.fna

# VCF of variant sites between BY and RM
VCFin=/home/mioverto/geneDrive/POS_files/RMxBY_ref_bcf.vcf
VCFnrm=${VCFin/.vcf/.nrm.vcf}
VCFgz=${VCFnrm/.nrm.vcf/.vcf.gz}

# Also create a chain file for correcting position indicies between references
chainRefToAlt=/home/mioverto/geneDrive/POS_files/RMvcf/BYtoRM.chain
# chainRMtoBY=/home/mioverto/geneDrive/POS_files/RMvcf/RMtoBY.chain

# Format VCF before creating chain file and masked sequence
bcftools norm -m +any -f $refIn $VCFin -Ov -o $VCFnrm

# bcftools requires compressed and indexed vcf file
bgzip -c $VCFnrm > $VCFgz
bcftools index $VCFgz

# Create alternative reference sequence and chain file
bcftools consensus -H A -f $refIn -c $chainRefToAlt $VCFgz > $newSeqOut


```
# If bam file lacks required metadata, this should help to troubleshoot
# Requires GATK

#Switch to Java 1.8 for GATK v4.X to function correctly. The long dir name may have to be updated when java is updated on the cluster
PATH=/usr/lib/jvm/java-1.8.0-openjdk-1.8.0.252.b09-2.el7_8.x86_64/bin:$PATH
# Path to GATK application
GATK=/home/mioverto/bin/gatk_4/gatk-package-4.1.8.0-local.jar

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
