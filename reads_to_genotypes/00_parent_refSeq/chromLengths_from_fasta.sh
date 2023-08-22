#!/bin/bash

while getopts 'R:o:' OPTION; do
   case "$OPTION" in
    R)
        refSeq="${OPTARG}"
        echo "Operating system $OPTARG"
        ;;

   o)
      outDir="$OPTARG"
      echo "Output directory $OPTARG"
      ;;


   ?)
        echo "script usage: $(basename \$0) [-v variantDirectory] [-l lineage] \
         [-f founderNameOrFileOfNames] [-o outputDirectory] [-R referenceSequence]"
    esac
done

# refSeq=${P1refSeq}
# outDir=${P1refSeq/_refseq.fna/_chrom_lengths.tsv}
cat ${refSeq} | \
awk '$0 ~ ">" {if (NR > 1) {print c;} c=0;printf substr($0,2,100) "\t"; } $0 !~ ">" {c+=length($0);} END { print c; }' > \
${outDir}
