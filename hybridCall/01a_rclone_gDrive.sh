#!/bin/bash
#PBS -A mioverto
#PBS -l nodes=1 
#PBS -l walltime=12:00:00

cd ${readsDir}

    # rclone was installed in home/mioverto/bin
    # Navigate to desired dir and use rclone copy

########
# Transfer read data from SKLAB_DATA Google Drive to scratch drive - run from TSCC in interactive

    #   MAseq1 - Pooled sequencing with Dutton lab 1/2 performed by SEG and MSO
rclone copy SKLAB_DATA:/2019_11_08_Dutton_Barcode_and_MutAccum_SG/MA_seq_data ${readsDir}

    #   MAseq2 - Pooled sequencing with Dutton lab 2/2 performed by SEG and MSO
rclone copy -P  --exclude Dutton_Data/** --exclude Reports/** --exclude Stats/** SKLAB_DATA:/2019_11_08_Kryazhimskiy_Barcode_and_MutAccum_SG ${readsDir}

    #   MAseq3 - Pooled sequencing with Alena 1/1 performed by MSO
rclone copy -P SKLAB_DATA:/2020_01_08_Kryazhimskiy_Robust_MutAcc_MO_AM/2020_01_08_MAseq_MO ${readsDir}


# rclone copy -P --include L1_* SKLAB_DATA:RAW_SEQ_DATA/2019_11_08_Dutton_Barcode_and_MutAccum_SG/MA_seq_data ./

