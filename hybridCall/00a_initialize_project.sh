#!/bin/bash

while getopts 's:p:b:a:c:t:O:T:' OPTION; do
    case "$OPTION" in
    s)
        opSys="$OPTARG"
        echo "System set to $OPTARG"
        ;;
    p)
        proj="$OPTARG"
        echo "Project name is $OPTARG"
        ;;
    b)
        baseDir="$OPTARG"
        echo "Base directory is set to $OPTARG"
        ;;
    a)
        appDir="$OPTARG"
        echo "Application directory is set to $OPTARG"
        ;;
    c)
        codeDir="$OPTARG"
        echo "Code directory is set to $OPTARG"
        ;;
    t)
        typeRef="$OPTARG"
        echo "The type reference strain is set to $OPTARG"
        ;;
    O)
        P1="$OPTARG"
        echo "Parent 1 is set to $OPTARG"
        ;;
    T)
        P2="$OPTARG"
        echo "Parent 2 is set to $OPTARG"
        ;;

    ?)
        echo "script usage: $(basename \$0) [-s computingSystem]
        [-p projectName] [-b baseDirectory] [-a applicationDir] [-c codeDir]
        [-t typeStrainName] [-O parentOneName] [-T parentTwoName]"
    esac
done



cd $baseDir
mkdir $proj
projDir=${baseDir}/${proj}
cd $projDir
sysParams=${projDir}/sysparams.txt
echo -e "${opSys}\n${proj}\n${baseDir}\n${appDir}\n${codeDir}" > $sysParams

# progDir=${codeDir}/hybridCall

### Still working on this part. Works but needs an if statement to prevent overwriting.

mkdir reads alignments metrics variants references metadata logs

cd ${projDir}/reads
mkdir sra_cache

cd ${projDir}/alignments
mkdir $P1 $P2

cd ${projDir}/metrics
mkdir reads alignments
cd reads
mkdir $P1 $P2
cd ../alignments
mkdir $P1 $P2

cd ${projDir}/variants
mkdir individuals lines final
cd individuals
mkdir $P1 $P2 
cd ../lines
mkdir $P1 $P2 
cd ../final
mkdir $P1 $P2 

cd ${projDir}/references
mkdir construction final POS_files
cd construction
mkdir reads alignments variants
cd reads
mkdir sra_cache $P1 $P2
cd ${projDir}/references/construction/alignments
mkdir $P1 $P2
cd ${projDir}/references/construction/variants
mkdir $P1 $P2
cd ${projDir}/references/final
mkdir $typeRef $P1 $P2

cd ${projDir}/logs
mkdir trim alignments variants

# binDir=/home/mioverto/bin

sysVars=${projDir}/sysvars.sh

rcloneDir=`find $appDir -type d -name "rclone*"`
echo "export rcloneApp=`find $rcloneDir -type f -name rclone`" >> ${sysVars}

sraDir=`find $appDir -type d -name "sratoolkit*"`
echo "export prefetch=${sraDir}/bin/prefetch" >> ${sysVars}
echo "export fasterq=${sraDir}/bin/fasterq-dump" >> ${sysVars}

cd $appDir/..
homeDir=`echo $(pwd)`
echo "export homeCache=`find ${homeDir} -type d -name "prefetc*"`/sra" >> ${sysVars}
echo "export fasterqCache=${projDir}/reads/sra_cache" >> ${sysVars}
echo "export pFasterqCache=${projDir}/references/construction/reads/sra_cache" >> ${sysVars}

trimDir=`find $appDir -type d -name "Trimmomatic*"`
echo "export trimApp=`find $trimDir -type f -name "trimmomatic*"`" >> ${sysVars}

FastqcDir=`find ${appDir} -maxdepth 2 -type d -name "FastQC*"`
echo "export FastQC=${FastqcDir}/fastqc" >> ${sysVars}

# echo "bedtools=${appDir}/bedtools" >> ${codeDir}/hybridCall/02_initialize_variables.sh

jvmDir=${appDir}/jvm
echo "export javaDir=`find $jvmDir -type d -name "java*"`/bin" >> ${sysVars}
# PATH=${javaDir}:$PATH

gatkDir=`find $appDir -maxdepth 1 -type d -name "gatk*"`
echo "export gatk=`find $gatkDir -type f -name "gatk*local.jar"`" >> ${sysVars}

# echo "export adapters=${metaDir}/NexteraPE-PE.fa" >> ${sysVars}