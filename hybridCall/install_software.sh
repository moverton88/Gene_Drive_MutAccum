#!/bin/bash

## rclone #####################################################################



## SRA  #######################################################################

cd ${HOME}/bin

wget --output-document sratoolkit.tar.gz https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-centos_linux64.tar.gz

tar -vxzf sratoolkit.tar.gz

export PATH=$PATH:${HOME}/bin/sratoolkit.3.1.0-centos_linux64/bin

which ${HOME}/bin/sratoolkit.3.1.0-centos_linux64/bin/fastq-dump

```
${HOME}/bin/sratoolkit.3.1.0-centos_linux64/bin/fastq-dump

# There are only a handfull of options that need to be enabled 
# to access public and controlled-access data in the cloud. 
# You will see a screen where you operate the buttons by pressing 
# the letter highlighted in red, or by pressing the tab-key until 
# the wanted button is reached and then pressing the space- or the enter-key.

## Enable the "Remote Access" option on the Main screen.

## In the "Cache" tab, enable "local file-caching"

## set the "Location of user-repository"
    # The repository directory needs to be set to an empty folder. 
    # This is the folder where prefetch will deposit the files.


# To start the configuration, run:


vdb-config -i

# Test that the software works
fastq-dump --stdout -X 2 SRR390728

## TRIMMOMATIC ################################################################

wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip
unzip Trimmomatic-0.39.zip
PATH=${HOME}/bin/Trimmomatic-0.39:$PATH


## htslib #####################################################################

cd ${HOME}/bin
wget https://github.com/samtools/htslib/releases/download/1.19.1/htslib-1.19.1.tar.bz2
tar xjvf htslib-1.19.1.tar.bz2
cd ${HOME}/bin/htslib-1.19.1    # and similarly for bcftools and htslib
./configure --prefix=${HOME}/bin --disable-bz2
make
make install
export PATH=${HOME}/bin/htslib-1.19.1:$PATH


## samtools ###################################################################

cd ${HOME}/bin
wget https://github.com/samtools/samtools/releases/download/1.19.2/samtools-1.19.2.tar.bz2
tar xjvf samtools-1.19.2.tar.bz2
cd ${HOME}/bin/samtools-1.19.2
./configure --prefix=${HOME}/bin --without-curses --disable-bz2
make
make install
export PATH=${HOME}/bin/samtools-1.19.2:$PATH  


## bcftools ###################################################################

cd ${HOME}/bin
wget https://github.com/samtools/bcftools/releases/download/1.19/bcftools-1.19.tar.bz2
tar xjvf bcftools-1.19.tar.bz2
cd ${HOME}/bin/bcftools-1.19
./configure --prefix=${HOME}/bin --without-curses --disable-bz2
make
make install
export PATH=${HOME}/bin/bcftools-1.19:$PATH  

## bedtools ###################################################################

cd ${HOME}/bin
# wget https://github.com/arq5x/bedtools2/releases/download/v2.31.1/bedtools-2.31.1.tar.gz
wget -r https://github.com/arq5x/bedtools2/releases/download/v2.28.0/bedtools
mv github.com/arq5x/bedtools2/releases/download/v2.28.0/bedtools ${HOME}/bin
# tar zxvf bedtools-2.31.1.tar.gz
# mv bedtools2 bedtools
# cd bedtools
# make LIBS='/PATH/TO/ZLIB/lib/libz.a'
export PATH=${HOME}/bin/:$PATH

## bwa ###################################################################

cd ${HOME}/bin
curl -L https://github.com/bwa-mem2/bwa-mem2/releases/download/v2.2.1/bwa-mem2-2.2.1_x64-linux.tar.bz2 \
  | tar jxf -
# wget https://sourceforge.net/projects/bio-bwa/files/bwa-0.7.16a.tar.bz2
# tar xjvf bwa-0.7.16a.tar.bz2
tar xjvf bwa-mem2-2.2.1_x64-linux.tar.bz2

## Java17 ###############################################################
cd ${HOME}/bin
# wget https://builds.openlogic.com/downloadJDK/openlogic-openjdk/8u402-b06/openlogic-openjdk-8u402-b06-linux-x64.tar.gz
wget https://builds.openlogic.com/downloadJDK/openlogic-openjdk/17.0.10+7/openlogic-openjdk-17.0.10+7-linux-x64.tar.gz
tar zxvf openlogic-openjdk-17.0.10+7-linux-x64.tar.gz
mv openlogic-openjdk-17.0.10+7-linux-x64 jdk17 
export PATH=${HOME}/bin/jdk17/bin:$PATH

## GATK #######################################################################
cd ${HOME}/bin
# wget https://github.com/broadinstitute/gatk/releases/download/4.5.0.0/gatk-4.5.0.0.zip
wget https://github.com/broadinstitute/gatk/releases/download/4.3.0.0/gatk-4.3.0.0.zip
unzip gatk-4.3.0.0.zip

## mummer 4 ###################################################################
cd ${HOME}/bin
wget https://github.com/mummer4/mummer/releases/download/v4.0.0rc1/mummer-4.0.0rc1.tar.gz
tar zxvf mummer-4.0.0rc1.tar.gz
cd mummer-4.0.0rc1
./configure
make
make install

## R ##########################################################################
cd ${HOME}/bin
wget https://cran.r-project.org/src/base/R-4/R-4.3.3.tar.gz
tar zxvf R-4.3.3.tar.gz
cd  R-4.3.3
module load cpu/0.17.3 gpu/0.17.3 gcc/10.2.0-mqbpsxf bzip2 curl/7.79.0-7bq2mgj
./configure --with-readline=no --with-x=no --enable-R-static-lib
make
make install

## fastqc ######################################################################
wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.12.1.zip
unzip fastqc_v0.12.1.zip
