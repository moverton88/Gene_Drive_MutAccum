

ssh mioverto@tscc-login.sdsc.edu


    # Start an interactive session
qsub -I -A svenkata

    # TSCC scratch drive to store bulk data 
cd /oasis/tscc/scratch/mioverto/mutAccum/reads

    # rclone was installed in home/mioverto/bin and authorized to be used everywhere
    # Navigate to desired dir and use rclone copy

########
# Transfer read data from SKLAB_DATA Google Drive to scratch drive - run from TSCC

    #   MAseq1 - Pooled sequencing with Dutton lab 1/2 performed by SEG and MSO
rclone copy SKLAB_DATA:/2019_11_08_Dutton_Barcode_and_MutAccum_SG/MA_seq_data ./MAseq1

    #   MAseq2 - Pooled sequencing with Dutton lab 2/2 performed by SEG and MSO
rclone copy -P  --exclude Dutton_Data/** --exclude Reports/** --exclude Stats/** SKLAB_DATA:/2019_11_08_Kryazhimskiy_Barcode_and_MutAccum_SG ./MAseq2

    #   MAseq3 - Pooled sequencing with Alena 1/1 performed by MSO
rclone copy -P SKLAB_DATA:/2020_01_08_Kryazhimskiy_Robust_MutAcc_MO_AM/2020_01_08_MAseq_MO ./MAseq3

