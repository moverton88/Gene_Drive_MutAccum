#!/bin/bash


module load cpu/0.17.3 
module load gcc/10.2.0-2ml3m2l
module load R

repeats_script=${codeDir}/repeats_from_masked_fasta.R
input_seq=${refSeqDir}/S288C/S288C_R64_refseq.fasta
output_bed=${refSeqDir}/S288C/repeats/S288C_repeats.bed
mkdir ${refDir}/final/S288C/repeats

Rscript $repeats_script -s $input_seq -b $output_bed -F TRUE
