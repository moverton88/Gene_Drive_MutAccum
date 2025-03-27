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

# refSeq=${refSeqDir}/${P_in}/${P_in}_refseq.fasta
# refSeq=/home/mioverto/geneDrive/refseq/RM/RM_refseq_UCSD_2020_v4.fna
# refSeq=${refSeqDir}/${P2}/${P2}_refseq.fasta
# lenOut=/home/mioverto/geneDrive/refseq/RM/RM_chrom_lengths.tsv
# lenOut=$(echo ${refSeq/_refseq.fasta/_chrom_lengths.tsv})
# lenOut="S288C_chrom_lengths.tsv"

cat ${refSeq} | awk '$0 ~ ">" {if (NR > 1) {print c;} c=0;printf substr($0,2,100) "\t"; } $0 !~ ">" {c+=length($0);} END { print c; }' | head -n 16 > ${lenOut}
