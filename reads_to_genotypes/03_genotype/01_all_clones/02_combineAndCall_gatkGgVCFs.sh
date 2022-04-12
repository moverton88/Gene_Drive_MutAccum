#!/bin/bash
#PBS -A mioverto
#PBS -l nodes=1 
#PBS -l walltime=5:00:00

# Combine VCF files of individual clones into one master VCF

PATH=/usr/lib/jvm/java-1.8.0-openjdk-1.8.0.252.b09-2.el7_8.x86_64/bin:$PATH
GATK=/home/mioverto/bin/gatk_4/gatk-package-4.1.8.0-local.jar

# lineage=all
if [ ${lineage} == all ] ; then
    export gVCFlist=${gVCFdir}/${lineage}_gVCFlist.list
    dir ${gVCFdir}/*${alignTag}.g.vcf > ${gVCFlist}

    export founderList=${gVCFdir}/${lineage}_founderList.list
    export founderList=${gVCFdir}/${lineage}_founderList.list
    for founder in ${gVCFdir}/*00*.g.vcf; do
        export index=$(basename "${founder}" .g.vcf)
        export id=${index:0:5}
        echo $id >> ${founderList}
    done
    else
        export gVCFlist=${gVCFdir}/${lineage}_gVCFlist.list
        dir ${gVCFdir}/${lineage}*${alignTag}.g.vcf > ${gVCFlist}
        export founderList=${gVCFdir}/${lineage}_founderList.list
        for founder in ${gVCFdir}/${lineage}*00*${alignTag}.g.vcf; do
            export index=$(basename "${founder}" .g.vcf)
            export id=${index:0:5}
            echo $id >> ${founderList}
        done
fi

# export gVCFlist=${gVCFdir}/anc_gVCFlist.list
# dir ${gVCFdir}/*00*.g.vcf >> ${gVCFlist}
# less ${gVCFlist}
# less ${founderList}

java -jar -Xmx16g $GATK CombineGVCFs  \
    -R $refSeq \
    --variant $gVCFlist \
    --founder-id ${founderList} \
    -O ${lineVCFdir}/${lineage}_multi.g.vcf.gz

if [ -z ${intervalFile} ]; then
    java -jar -Xmx16g $GATK GenotypeGVCFs \
        -R $refSeq \
        -V ${lineVCFdir}/${lineage}_multi.g.vcf.gz \
        --founder-id ${founderList} \
        -O ${VCFout}
fi

if [ ! -z ${intervalFile} ]; then
    java -jar -Xmx16g $GATK GenotypeGVCFs \
        -R $refSeq \
        -V ${lineVCFdir}/${lineage}_multi.g.vcf.gz \
        --founder-id ${lineage}00 \
        -L $intervalFile \
        -all-sites true \
        -O ${VCFout}
fi

```
export founderList=${gVCFdir}/A_founderList.list

for founder in ${gVCFdir}/${lineage}*00*.g.vcf; do
    export index=$(basename "${founder}" .g.vcf)
    export id=${index:0:5}
    echo $id >> ${founderList}
done


export gVCFlist=${gVCFdir}/anc_gVCFlist.list
dir ${gVCFdir}/*00*.g.vcf >> ${gVCFlist}
less ${gVCFlist}

java -jar -Xmx16g $GATK CombineGVCFs  \
    -R $refseq \
    --variant $gVCFlist \
    -O ${lineVCFdir}/anc_multi.g.vcf.gz

VCFout=${finalVCFdir}/A_RMaRMc_allVar.vcf

java -jar -Xmx32g $GATK GenotypeGVCFs \
    -R $refseq \
    --founder-id ${founderList} \
    -V ${lineVCFdir}/A_multi.g.vcf.gz \
    -O ${VCFout}





alignRef=RM
callRef=RM
lineage=N_A
# refseq=/home/mioverto/geneDrive/refseq/BY/S288C_R64_refseq.fna
refseq=/home/mioverto/geneDrive/refseq/RM/RM_refseq_UCSD_2020_v4.fna
lineVCFdir=/oasis/tscc/scratch/mioverto/geneDrive/dualRef/${alignRef}_aligned/variants/gVCFs/${callRef}_call/lineageVCFs
finalVCFdir=/oasis/tscc/scratch/mioverto/geneDrive/dualRef/${alignRef}_aligned/variants/LOHvcfs/${callRef}_call
# lineage=N_B
VCFout=${finalVCFdir}/${lineage}_${alignRef}align_${callRef}call.test.vcf
intervalFile=/home/mioverto/geneDrive/POS_files/${callRef}_POS.bed


java -jar $GATK GenotypeGVCFs \
   -R $refseq \
   -V ${lineVCFdir}/${lineage}_multi.g.vcf.gz \
   -L $intervalFile \
   -all-sites true \
   -O ${VCFout}



```