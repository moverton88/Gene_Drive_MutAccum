#!/bin/bash

###############################################################################
# Run in interactive mode #
###############################################################################

# Takes a planfile of old and new filenames and renames files in rename_dir 
# to the new file names old_name must be a filename without the extention. 
#It must identify a single sample, but will rename paired-read files

export rename_dir=/oasis/tscc/scratch/mioverto/geneDrive/reads/raw/
export plan_in=/home/mioverto/code/rename/masterPlanFile.csv
export plan_file=/home/mioverto/code/rename/masterPlanFile_crct.csv
# export seqID=1
cd $rename_dir

# find ./ -name "L218""_*R2*"

# unset old_name=L218
# unset new_name=F_B14

# remove carriage return from masterPlanFile.csv > masterPlanFile_v2.csv 
# must use control-v then control-m to enter ^M
sed -e “s/^M//” $plan_in > $plan_file
# in vim
:%s/^M//g

IFS=,
[ ! -f "$plan_file" ] && { echo "$plan_file file not found"; exit 99; }
while read old_name new_name
do
    #echo "searching for $old_name to replace with $new_name"
#done < "$plan_file"
    oldR1=$(find ./ -name "${old_name}*_R1*")
    oldR2=$(find ./ -name "${old_name}*_R2*")
    if [ -f "$oldR1" ]; then
        echo "$oldR1 found"
        base=$(basename ${oldR1})
        sfx="${base#*.}"
        newR1=$(echo "${new_name}_R1.${sfx}")
        echo "will rename ${oldR1} to ${newR1}"
        mv $oldR1 $newR1
    fi
# done < "$plan_file"
    if [ -f "$oldR2" ]; then
        newR2=$(echo ${new_name}_R2.${sfx})
        echo "will rename ${oldR2} to ${newR2}"
        mv ${oldR2} ${newR2}
    fi
done < "$plan_file"

```
# manual correction

old_nm=L24_S32_L002_R2_001.fastq.gz
new_nm=N_B12_R2.fastq.gz
mv $old_nm $new_nm

```