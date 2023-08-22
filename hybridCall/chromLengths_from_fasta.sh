#!/bin/bash

while getopts 'R:o:' OPTION; do
   case "$OPTION" in
    R)
        refSeq="${OPTARG}"
        echo "Operating system $OPTARG"
        ;;

   o)
      lenOut="$OPTARG"
      if [[ ${lenOut} == "" ]]; then
         echo "no length file given, creating from fasta"
         lenOut=${refSeq/_refseq.fna/_chrom_lengths.tsv}

      else
         echo "length file is $lenOut"
      fi
      ;;


   ?)
        echo "script usage: $(basename \$0) [-R referenceSequence] [-o chromLengthFileName]"
    esac
done

# refSeq=${P1refSeq}
# refSeq=/home/mioverto/geneDrive/refseq/RM/RM_refseq_UCSD_2020_v4.fna
# refSeq=/oasis/tscc/scratch/mioverto/LOH_methods/Overton_etal/references/final/RM/RM_refseq.fasta
# lenOut=/home/mioverto/geneDrive/refseq/RM/RM_chrom_lengths.tsv
# lenOut=${refSeq/_refseq.fasta/_chrom_lengths.tsv}
# lenOut=${refSeq/_refseq${tag}.fna/_chrom_lengths${tag}.tsv}

cat ${refSeq} | \
awk '$0 ~ ">" {if (NR > 1) {print c;} c=0;printf substr($0,2,100) "\t"; } $0 !~ ">" {c+=length($0);} END { print c; }' > \
${lenOut}
