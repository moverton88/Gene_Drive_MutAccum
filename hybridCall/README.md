---
title: "hybridCall: Dual reference alignment, variant calling, and genotype reconciliation"
author: "Michael Overton"
output: html_notebook
---

<!--########################################################################-->
# hybridCall
<br>
hybridCall is a shell-based set of scripts that constructs a pair of parental 
reference sequences, and then aligns and calls variants against each reference 
sequence for a set of samples.
<br>
The programatic structure of hybridCall relies on three main directories: 
project, code, and applications. The project directory is structured as follows:
<br>
(P1 = Parent 1 strain name, P2 = Parent 2 strain name)
<br> <br>
projDir<br>
|-references<br>
 |-construction<br>
  |-reads<br>
   |-P1<br>
   |-P2<br>
  |-alignments<br>
  |-variants<br>
 |-final<br>
|- reads<br>
|
|
|
|

<br><br>
# 00 Create the project and directory structure

The first step in executing this pipeline is to create the requred directories and meta-files. To do this, 
the script "00a_initialize_project.sh" is run with the following options:
<br>
The script will create a project directory in the base directory, create a "sysparams.txt" file to 
store directory information, construct the remaining directory structure, and create a "sysvars.sh"
that stores the locations of the applications required for the pipeline. Before running this script, 
the following steps must be completed:
<br>
Install the required applications in the application directory:
- rclone
- SRA toolkit
- Trimmomatic
- Fast QC
- Java JVM v1.8
- GATK v4.0+
- mummer v4.0+

<br>
Place the Hybrid Call folder in the code directory.
<br><br>
Once these steps are complete, run 00a_initialize_project.sh with the following options.
    -s operating system ("tscc" for TSCC remote cluster or "other" for any other OS)
    -p project name
    -b base directory (directory in which to create project)
    -a application directory (where the required applications, such as bwa and GATK, are located)
    -c code directory (where the Hybrid Call scripts are located)
    -t Type reference strain name
    -O Parent 1 strain name
    -T Parent 2 strain name
<br><br>
# 00b Install required software packages

<br><br>
## rclone

<br><br>
## SRA toolkit

<br><br>
## Trimmomatic

<br><br>
## Fast Qc

<br><br>
## Java JVM

<br><br>
## GATK

<br><br>
## mummer

<br><br>

# 01 Obtain sequencing reads/contigs for each parent and all samples

<p>
Before processing the parental and sample sequences, they need to be placed 
into their respective directories. The parental sequences need to be placed 
into the {project}/references/construction/reads/{parent name}/ directory and 
the sample reads into the {project}/reads/ directory. You can either manually 
place them into the appropriate folder or download them from the NCBI SRA 
repository using SRA toolkit. 
<p\>
<p>
To use SRA toolkit, place a .txt file containing the NCBI Accession IDs separated 
by new lines into the {project}/metadata directory. Then, run 
"01_download_SRA_reads.sh" with the paramters:<br>
-A file with accession IDs<br>
-C SRA cache directory (*./references/construction/reads/sra_cache/ -or- #./reads/sra_cache/)<br>
-O output directory (*./references/construction/reads/{parent name}/ -or- #./reads/ -or- ^./reads/raw)<br>
<blockquote> 
    * Parent directories<br>
    # Sample directories<br>
    ^ If reads will require trimming<br>
</blockquote>
<p>
<br>

# 02 Read pre-processing

Read file names are often long and unintuitive. Therefore, a renaming script, "02a_rename_fastq.sh" 
is included. This script takes in a comma-separated value file with each line containing an 
[oldName, newName] pair, and a directory in which the read files are located. The script will 
search for pairs of files that match the oldName and rename the files {newName}_R1.fastq and 
{newName}_R2.fastq.

Often, the ends of sequencing reads contain adapter sequences and/or low-quality base calls and must 
be trimmed. This is achieved with the "02b_trim_reads.sh" script, which runs the reads through 
Trimmomatic. The script takes in the input reads file (either a single-end read file or pair of 
paired-end read files), an option for single-end (SE) or paired-end (PE), the output directory,
and all of the options for Trimmomatic (see manual).

# 03 Construct parental reference sequences

Place the Type strain reference sequence must be placed into the {project}/references/final/{Type name}/
directory.

Once the parent reads/contigs and Type reference sequence are placed in their directories, run either
"03_construct_refSeq.sh" or "03_construct_refSeq_contigs.sh" depending on the type of input sequences. 


# 04 Sample read alignment and variant calling

Once the parental reference sequences have been constructed and the sample reads are prepared,
alignment and variant calling can be carried out. This is done by running the "00_run_indv_pipe.sh".
This script takes in the following information:

- Project name 
- Sample name
- Type strain name
- Parent name
- Variables file



