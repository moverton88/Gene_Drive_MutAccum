#!/bin/bash

###############################################################################
# Run in interactive mode #
###############################################################################

# Takes a planfile of old and new filenames and renames files in rename_dir 
# to the new file names old_name must be a filename without the extention. 
# It must identify a single sample, but will rename paired-read files

export proj=Dutta_etal_2021
export rename_dir=/oasis/tscc/scratch/mioverto/LOH_methods/${proj}/parents/reads/raw
# export rename_dir=/oasis/tscc/scratch/mioverto/LOH_methods/Pankajam_etal_2020/reads
# export rename_dir=/oasis/tscc/scratch/mioverto/geneDrive/reads/raw/

export plan_dir=/home/mioverto/LOH_methods/rename
# export plan_dir=/home/mioverto/code/rename

export plan_name=${proj}_planfile.csv
# export plan_name=Pankajam_etal_2020_planfile.csv
# export plan_name=masterPlanFile.csv

export plan_file=${plan_dir}/${plan_name}
export plan_crct=${plan_file/.csv/_crct.csv}

cd $rename_dir

# remove carriage return from masterPlanFile.csv > masterPlanFile_v2.csv 
# must use control-v then control-m to enter ^M
# sed -e “s/^M//” $plan_file > $plan_crct


# in vim
:%s/^M//g

# oldR1=SRR11460413_1.fastq
# old_name=SRR11460413
# new_name=L_1

while IFS="," read -r accNum ID || [[ -n "$accNum" ]]; do
#    echo "searching for $accNum"
#    echo "to replace with $ID"
# done < "$plan_file"
    acc=$(echo "$accNum")
    id=$(echo "$ID")
    oldR1=$(find ./ -name "${acc}*1*")
    oldR2=$(find ./ -name "${acc}*2*")
    if [ -f "$oldR1" ]; then
        base=$(basename ${oldR1})
        sfx="${base#*.}"
        newR1=$(echo "${id}_R1.${sfx}")
        echo "$oldR1 found, will rename to ${newR1}"
        mv $oldR1 $newR1
    fi
# done < "$plan_file"
    if [ -f "$oldR2" ]; then
        newR2=$(echo ${id}_R2.${sfx})
        echo "$oldR2 found, will rename to ${newR2}"
        mv ${oldR2} ${newR2}
    fi
# done < "$plan_crct"
done < "$plan_file"

```
# manual correction
acc=SRR12442664
find ./ -name "${acc}*"

old_nm=L24_S32_L002_R2_001.fastq.gz
new_nm=N_B12_R2.fastq.gz
mv $old_nm $new_nm

```
SRR12442617
SRR12442618
SRR12442626
SRR12442637
SRR12442648
SRR12442667
SRR12442671
SRR12442672
SRR12442679
SRR12442680
SRR12442686
SRR12442688
SRR12442690
SRR12442691
SRR12442696
SRR12442697
SRR12442698
SRR12442699
SRR12442702
SRR12442703
SRR12442705

