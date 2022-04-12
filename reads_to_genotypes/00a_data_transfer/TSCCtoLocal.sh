

# Define from and to directories

## Bam file of RM scaffolds aligned to BY reference
fromDir=/oasis/tscc/scratch/mioverto/geneDrive/refseq/RM/RM_refseq.crct.bam.bai
toDir=./

## Final RM reference sequence
fromDir=/home/mioverto/geneDrive/refseq/RM/RM_refseq_UCSD_2020_v4.fna

## Composite RM ref RMxBY vcf
fromDir=/home/mioverto/geneDrive/POS_files/RMxBY_comp_ref_bcf.vcf
fromDir=/home/mioverto/geneDrive/POS_files/RMxBY_comp_ref_HC.vcf

## Bam file of clone reads mapped to either reference
fromDir=/oasis/tscc/scratch/mioverto/geneDrive/dualRef/RM_aligned/bam/DeDup/
fromDir=/oasis/tscc/scratch/mioverto/geneDrive/dualRef/BY_aligned/bam/DeDup/
fromDir=/oasis/tscc/scratch/mioverto/geneDrive/dualRef/BY_aligned/variants/gVCFs/BY_call/realign_bam/
fromDir=/oasis/tscc/scratch/mioverto/geneDrive/dualRef/RM_aligned/variants/gVCFs/RM_call/realign_bam/

fromDir=/oasis/tscc/scratch/mioverto/geneDrive/Cas9/bam/DeDup

inc="F_D*.ba*"

## Final VCFs of lineages aligned to BYm and called against the RM reference
fromDir=/oasis/tscc/scratch/mioverto/geneDrive/ambiRef/variants/LOHvcfs/RM_call
toDir=./

## VCF of ancestors
fromDir=/oasis/tscc/scratch/mioverto/geneDrive/dualRef/BY_aligned/variants/allVarVcfs
toDir=./
inc="anc_pool*.g.vcf"


## Final VCFs of lineages aligned to the BY reference and called against the BY reference
fromDir=/oasis/tscc/scratch/mioverto/geneDrive/dualRef/BY_aligned/variants/allVarVcfs/filtered
## Final VCFs of lineages aligned to the RM reference and called against the BY reference
fromDir=/oasis/tscc/scratch/mioverto/geneDrive/dualRef/RM_aligned/variants/allVarVcfs/filtered
## Final VCFs of lineages aligned to the Cas9 gRNA reference
fromDir=/oasis/tscc/scratch/mioverto/geneDrive/Cas9/variants/allVarVcfs/

toDir=./
inc="*filtered.vcf"
inc="all*default_allVar.vcf"

toDir=./
inc="*.vcf"

## Final table of lineages aligned to the BY reference and called against the BY reference
fromDir=/oasis/tscc/scratch/mioverto/geneDrive/dualRef/RM_aligned/variants/allVarVcfs/tables
inc="*.tsv"

## VCF from bcftools call
fromDir=/oasis/tscc/scratch/mioverto/geneDrive/dualRef/BY_aligned/variants/bcftools
inc="N_E*_bcf_local.vcf"

## RMxBY vcf file
fromDir=/home/mioverto/geneDrive/POS_files/RMvcf/RMxBY_ref_noMit.vcf

## Chain file for translating positional indicies
fromDir=/home/mioverto/geneDrive/POS_files/RMvcf/RMxBY_ref_bcf.chain

## No coverage files
alignRef=RM
fromDir=/oasis/tscc/scratch/mioverto/geneDrive/dualRef/${alignRef}_aligned/variants/gVCFs/indvSingle/noCover
toDir=./
inc="*.txt"

fromDir=/oasis/tscc/scratch/mioverto/geneDrive/dualRef/BY_aligned/variants/allVarVcfs/vqsr
toDir=./
inc="*.R"

# Picard read depth stats
fromDir=/oasis/tscc/scratch/mioverto/geneDrive/dualRef/BY_aligned/bam/DeDup
toDir=./
inc="*.txt"

# Picard read duplicates stats
fromDir=/oasis/tscc/scratch/mioverto/geneDrive/dualRef/BY_aligned/bam/picard_metrics
toDir=./
inc="*.txt"

# Depth stats
fromDir=/oasis/tscc/scratch/mioverto/geneDrive/dualRef/BY_aligned/bam/depth_metrics
fromDir=/oasis/tscc/scratch/mioverto/geneDrive/dualRef/RM_aligned/bam/depth_metrics

toDir=./
inc="*.txt"


# Retrieve files using rclone
rclone copy TSCC:$fromDir --include $inc $toDir
rclone copy TSCC:$fromDir $toDir

rclone config
