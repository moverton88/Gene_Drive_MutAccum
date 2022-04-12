########
# REFERENCE SEQUENCES
# scp ~/SK_Lab/PhD_Projects/geneDrive/refseq/BY_ref/S288C_R63_RM.fasta mioverto@tscc-login.sdsc.edu:/oasis/tscc/scratch/mioverto/data/refseq/RM_ref
scp ~/SK_Lab/PhD_Projects/geneDrive/refseq/RM_ref/RM11-1a_UCI_2019.fna mioverto@tscc-login.sdsc.edu:/home/mioverto/geneDrive/refseq/RM
scp ~/SK_Lab/PhD_Projects/geneDrive/refseq/BY_ref/S288C_R64_refseq.fna mioverto@tscc-login.sdsc.edu:/home/mioverto/geneDrive/refseq/BY

    # For making RM reference sequence
scp ~/SK_Lab/PhD_Projects/geneDrive/refseq/RM_ref/UCI_scaffolds/RM11-1a_UCI_2019.fna \
mioverto@tscc-login.sdsc.edu:/home/mioverto/geneDrive/refseq/RM/UCI_scaffolds

scp ~/SK_Lab/PhD_Projects/geneDrive/refseq/RM_ref/other_refs/*.fna \
mioverto@tscc-login.sdsc.edu:/home/mioverto/geneDrive/refseq/RM/other_refs

scp ~/SK_Lab/PhD_Projects/geneDrive/refseq/POS_files/RMxBY_ref_bcf.vcf \
mioverto@tscc-login.sdsc.edu:/home/mioverto/geneDrive/POS_files

    # Final RM reference sequence
scp ~/SK_Lab/PhD_Projects/geneDrive/refseq/RM_ref/RM_refseq_UCSD_2020_v4.fna mioverto@tscc-login.sdsc.edu:/home/mioverto/geneDrive/refseq/RM

    # Final BYm reference sequence
scp ~/SK_Lab/PhD_Projects/geneDrive/refseq/RM_ref/S288C_R64_RMmasked.fna mioverto@tscc-login.sdsc.edu:/home/mioverto/geneDrive/refseq/BYm

    # Cas9 gRNA constructs
# scp ~/SK_Lab/PhD_Projects/geneDrive/refseq/Cas9_gRNA_construct.fasta mioverto@tscc-login.sdsc.edu:/home/mioverto/geneDrive/refseq/
scp -r ~/SK_Lab/PhD_Projects/geneDrive/refseq/Drive/ mioverto@tscc-login.sdsc.edu:/home/mioverto/geneDrive/refseq/Drive/

########
# SCRIPT FILES - from local disc to TSCC - run from local drive

scp ~/Downloads/rclone-v1.57.0-linux-amd64.zip mioverto@tscc-login.sdsc.edu:/home/mioverto/bin

    ## - Rename read files ###################################################
# scp ~/SK_Lab/PhD_Projects/geneDrive/code/Rename/rename_fastq_bash_012720.sh mioverto@tscc-login.sdsc.edu:/home/mioverto/code/rename
# scp ~/SK_Lab/PhD_Projects/geneDrive/code/Rename/rename_fastq_030620.py mioverto@tscc-login.sdsc.edu:/home/mioverto/code/rename
scp ~/SK_Lab/PhD_Projects/geneDrive/code/rename/renameFastqFiles.sh mioverto@tscc-login.sdsc.edu:/home/mioverto/code/rename

    ### planfile
scp ~/SK_Lab/PhD_Projects/geneDrive/code/rename/masterPlanFile.csv mioverto@tscc-login.sdsc.edu:/home/mioverto/code/rename

    ## - Trim reads
scp ~/SK_Lab/PhD_Projects/geneDrive/code/trim/submitTrimReads.sh mioverto@tscc-login.sdsc.edu:/home/mioverto/code/trim
scp ~/SK_Lab/PhD_Projects/geneDrive/code/trim/trimReads_V2.sh mioverto@tscc-login.sdsc.edu:/home/mioverto/code/trim

    ## - Bowtie align and output bam
# scp ~/SK_Lab/PhD_Projects/geneDrive/code/align/submit_alignToBam_v1.sh mioverto@tscc-login.sdsc.edu:/home/mioverto/code/align
# scp ~/SK_Lab/PhD_Projects/geneDrive/code/align/alignToBam_v1.sh mioverto@tscc-login.sdsc.edu:/home/mioverto/code/align
# scp ~/SK_Lab/PhD_Projects/geneDrive/code/align/alignToBam_local.sh mioverto@tscc-login.sdsc.edu:/home/mioverto/code/align
scp ~/SK_Lab/PhD_Projects/geneDrive/code/reads_to_genotypes/02_reads_align/alignAndMarkDups.sh mioverto@tscc-login.sdsc.edu:/home/mioverto/code/align

    ## - DeDuplicate reads separately if needed
scp ~/SK_Lab/PhD_Projects/geneDrive/code/align/deDup_reads.sh mioverto@tscc-login.sdsc.edu:/home/mioverto/code/align

    ## -- Apply BSQR to bam files
 scp ~/SK_Lab/PhD_Projects/geneDrive/code/align/BQSR_calc.sh mioverto@tscc-login.sdsc.edu:/home/mioverto/code/align   

    ## -- Get depth stats
 scp ~/SK_Lab/PhD_Projects/geneDrive/code_v2/reads_to_genotypes/02_reads_align/calcDepthStats.sh mioverto@tscc-login.sdsc.edu:/home/mioverto/code/align   

    ## - HaplotypeCaller call all genotypes for pooled ancestor
scp -r ~/SK_Lab/PhD_Projects/geneDrive/code_v2/reads_to_genotypes/03_genotype/00b_refSeqAndPOS/call_AncHC.sh \
    mioverto@tscc-login.sdsc.edu:/home/mioverto/code/variants

scp -r ~/SK_Lab/PhD_Projects/geneDrive/code_v2/variants/ancestor_calling/callAncGVCF.sh mioverto@tscc-login.sdsc.edu:/home/mioverto/code/variants

    ## - HaplotypeCaller gVCF call full genome VCFs for each individual, combine into multisample gVCF, and call final variants.
# scp -r ~/SK_Lab/PhD_Projects/geneDrive/code/variants/callGvcfPipe_v3.sh mioverto@tscc-login.sdsc.edu:/home/mioverto/code/variants
scp -r ~/SK_Lab/PhD_Projects/geneDrive/code_v2/reads_to_genotypes/03_genotype/01_all_clones/01_call_gVCF_gatkHC.sh \
    mioverto@tscc-login.sdsc.edu:/home/mioverto/code/variants

scp -r ~/SK_Lab/PhD_Projects/geneDrive/code_v2/reads_to_genotypes/03_genotype/01_all_clones/01_call_gVCF_gatkHC_default.sh \
    mioverto@tscc-login.sdsc.edu:/home/mioverto/code/variants


    ## - Combine gVCFs and GenotypegVCFs pipeline - combine individual gVCFS into multisample gVCF, and call final variants.
scp -r ~/SK_Lab/PhD_Projects/geneDrive/code/variants/lineageCallToVcf.sh mioverto@tscc-login.sdsc.edu:/home/mioverto/code/variants
scp -r ~/SK_Lab/PhD_Projects/geneDrive/code/variants/submit_lineageCallPipe.sh mioverto@tscc-login.sdsc.edu:/home/mioverto/code/variants

scp -r ~/SK_Lab/PhD_Projects/geneDrive/code/variants/lineageCallToVcf_allVar.sh mioverto@tscc-login.sdsc.edu:/home/mioverto/code/variants
scp -r ~/SK_Lab/PhD_Projects/geneDrive/code/variants/lineageCallToVcf_allVarXL.sh mioverto@tscc-login.sdsc.edu:/home/mioverto/code/variants

scp -r ~/SK_Lab/PhD_Projects/geneDrive/code/reads_to_genotypes/03_genotype/01_all_clones/02_combineAndCall_gatkGgVCFs.sh \
    mioverto@tscc-login.sdsc.edu:/home/mioverto/code/variants

    # Interval files for calling against BY- and RM-reference sequences
scp -r ~/SK_Lab/PhD_Projects/geneDrive/refseq/POS_files/BY_call.intervals \
    mioverto@tscc-login.sdsc.edu:/home/mioverto/geneDrive/POS_files

scp -r ~/SK_Lab/PhD_Projects/geneDrive/refseq/POS_files/BY_POS.bed \
    mioverto@tscc-login.sdsc.edu:/home/mioverto/geneDrive/POS_files

scp -r ~/SK_Lab/PhD_Projects/geneDrive/refseq/POS_files/RM_POS.bed \
    mioverto@tscc-login.sdsc.edu:/home/mioverto/geneDrive/POS_files

    # RMxBY position file for filtering
scp -r ~/SK_Lab/PhD_Projects/geneDrive/refseq/POS_files/RMxBY_variants.vcf \
    mioverto@tscc-login.sdsc.edu:/home/mioverto/geneDrive/POS_files

scp -r ~/SK_Lab/PhD_Projects/geneDrive/refseq/POS_files/RMxBY_ref_bcf_noMit.vcf \
    mioverto@tscc-login.sdsc.edu:/home/mioverto/geneDrive/POS_files

# - Call genotypes with bcftools
scp -r ~/SK_Lab/PhD_Projects/geneDrive/code/variants/bcftools_call_VCF_020821.sh mioverto@tscc-login.sdsc.edu:/home/mioverto/code/variants

    # Bed file with repeat and RMxBY positions
scp -r ~/SK_Lab/PhD_Projects/geneDrive/refseq/POS_files/rpt_RMxBY_toMask.bed \
    mioverto@tscc-login.sdsc.edu:/oasis/tscc/scratch/mioverto/geneDrive/refseq/POS_files

# - Call deNovo mutations with GATK HC and bed file containing repeat and RMxBY sites to mask
scp -r ~/SK_Lab/PhD_Projects/geneDrive/code/variants/deNovo/submit_deNovoCall.sh mioverto@tscc-login.sdsc.edu:/home/mioverto/code/variants
scp -r ~/SK_Lab/PhD_Projects/geneDrive/code/variants/deNovo/deNovoCallVcf.sh mioverto@tscc-login.sdsc.edu:/home/mioverto/code/variants
scp -r ~/SK_Lab/PhD_Projects/geneDrive/code/variants/deNovo/Chrom.list mioverto@tscc-login.sdsc.edu:/home/mioverto/code/variants

scp -r ~/SK_Lab/PhD_Projects/geneDrive/code/variants/deNovo/deNovoCombineCall.sh mioverto@tscc-login.sdsc.edu:/home/mioverto/code/variants

    # - Full align-to-VCF pipeline
scp -r ~/SK_Lab/PhD_Projects/geneDrive/code/alignToVcfPipe/Submit_MO-pipe.sh mioverto@tscc-login.sdsc.edu:/home/mioverto/code/fullPipe
scp -r ~/SK_Lab/PhD_Projects/geneDrive/code/alignToVcfPipe/alignToVcfPipe_v4.sh mioverto@tscc-login.sdsc.edu:/home/mioverto/code/fullPipe
scp -r ~/SK_Lab/PhD_Projects/geneDrive/code/alignToVcfPipe/callToVcfPipe.sh mioverto@tscc-login.sdsc.edu:/home/mioverto/code/fullPipe
