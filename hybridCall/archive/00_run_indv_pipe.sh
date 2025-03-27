#!/bin/bash
#PBS -A mioverto
#PBS -l nodes=1 
#PBS -l walltime=48:00:00

###############################################################################
# 00_run_indv_pipe.sh
# Sequentially runs the scripts to align sample reads to a chosen reference 
# sequence and call variants.
#*****************************************************************************#

while getopts 'v:s:' OPTION; do
   case "$OPTION" in
   v)
      variablesFile="$OPTARG"
      echo "The file containing all variables is $OPTARG"
      ;;

   s)
      opSys="${OPTARG}"
      echo "Operating system $OPTARG"
      ;;

   ?)
        echo "script usage: $(basename \$0) [-R1 readsfile1] [-s operatingSystem]"
    esac
done

if [[ ${opSys} == "tscc" ]]; then
   echo "loading packages from modules"
      module load bwa
      module load samtools    
      module load bcftools
      module load bedtools
   else
   echo "loading packages from files"
      bwa=${bwaApp}
      samtools=${samtoolsApp}
      bcftools=${bcftoolsApp}
      bedtools=${bedtoolsApp}
fi
