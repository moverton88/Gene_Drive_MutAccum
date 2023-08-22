---
title: "Hybrid Call: Dual reference alignment and variant calling"
author: "Michael Overton"
output: html_notebook
---


Hybrid call is a shell-based set of scripts that constructs a pair of parental reference sequences, 
and then aligns and calls variants against each reference sequence for a set of samples. 


# 00 Create the project and directory structure

The first step in this pipeline is to create the requred directories and meta-files. To do this, 
the script "00a_initialize_project.sh" is run with options for the operating system (whether
the TSCC remote cluster or another system is being used), project name, base directory 
(directory in which to create project), application directory (where the required applications, 
such as bwa and GATK, are located), the code directory (where the Hybrid Call scripts are located),
the type strain name, and the name of the two parental strains. The script will create a project
directory in the base directory, create a "sysparams.txt" file to store directory information,
construct the remaining directory structure, and create a "sysvars.sh" that stores the locations
of the applications required for the pipeline. Before running this script, the following steps 
must be completed:

Install the required applications in the application directory:
rclone
SRA toolkit
Trimmomatic
Fast QC
Java JVM v1.8
GATK v4.0+
mummer v4.0+

Place the Hybrid Call folder in the code directory

Once these steps are complete, run 00a_initialize_project.sh with the appropriate options.

# 01 Obtain sequencing reads/contigs for each parent and all samples

Before processing the parental and sample sequences, they need to be placed into their
respective directories. The parental sequences need to be placed in the 
{project}/references/construction/reads/{parent name}/ directory and the sample reads into
the --- directory. You can either manually place them into the appropriate folder or 
download them from the NCBI SRA repository using SRA toolkit. 

For using the SRA toolkit, place a .txt file with the NCBI Accession ID into the {project}/metadata directory,
then run 01_download_SRA_reads.sh with the path/name of the Accession list, the sra cache directory
(for parental sequences: {project}/references/construction/reads/sra_cache/ and sample sequences: {project}/reads/sra_cache), and the output directory 
(parent:{project}/references/construction/reads/{parent name}/ and sample: {project}/reads/ -or- {project}/reads/raw if read trimming will be applied).

# 02 Construct parental reference sequences

Place the Type strain reference sequence must be placed into the {project}/references/final/{Type name}/
directory.

Once the parent reads/contigs and Type reference sequence are placed in their directories, run either
"03_construct_refSeq.sh" or "03_construct_refSeq_contigs.sh" depending on the type of input sequences. 

# 02 Align and 

Before cons
