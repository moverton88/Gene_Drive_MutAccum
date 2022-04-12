# Search gVCF for hom ref blocks with 0 coverage and count number of 
# positions in these blocks

# Run in interactive mode


alignRef=RM
export gVCFdir=/oasis/tscc/scratch/mioverto/geneDrive/dualRef/${alignRef}_aligned/variants/gVCFs/indvSingle

for gVCF in ${gVCFdir}/*.vcf; do
    gVCFname=$(basename "${gVCF}" .g.vcf)
    index=${gVCFname:0:8}
    echo ${index}
    grep 0/0:0 ${gVCF} | grep -v NC_001224.1 | grep -v DP= > ${gVCFdir}/noCover/${index}_noCover.txt
done
