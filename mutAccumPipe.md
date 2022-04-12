# Introduction  
   
## Project Background
Endonucleases, such as Cas-9, interact directly with DNA in target cells, often so that precise breaks can be introduced as a first step in altering sequence. Both the binding and cutting capacities of endonucleases have the potential to increase the rate of mutations not only at the target site, but throughout the genome. The change in mutation rates for Cas-9 transfected cells has not been well characterized, a critical step in improving the design of gene editing tools, such as gene drives. We introduced a Cas-9 based gene drive into yeast cells, and evolved them along mutation accumulation lineages to determine whether the presence of the entire drive and the Cas-9 protein alone can increase the rate of break repair recombination (gene conversion) and de novo mutations.

## Experiment Outline
Clones from a cross between laboratory (BY) and vineyard (RM) yeast strains were tranformed with one of three treatments: a negative control containing empty constructs, a half-drive consisting of a Cas-9 plasmid and integrated empty guide RNA construct, and a full-drive that includes the Cas-9 plasmid and integrated guide RNA targeting the ADE2 locus. Successful construct clones were used to seed mutation accumulation lines as follows: 8 individual transformants were recovered for each treatment and used to seed a founding population, and isolated clones from each population seeded 10-12 replicate mutation accumulation lines. Each replicate was diluted and spread onto a YPD plate, allowed to grow 48h, and a single colony was diluted and spread on the next plate. This was repeated for ~880 generations. Final colonies were grown in liquid YPD and sequenced to ~50x coverage. 

> *Sample Names*  
The naming convention in the analysis was changed from the original experiment to specify treatment (N=no drive, H=half-drive, F=full-drive), founder group (A-H), and replicate lineage (01-13) information. The new names appear as N_A00. N for No Drive, A for the first lineage, and 01 for the final evolved clone in the first replicate. A 00 in the replicate position indicates the ancestral clone of that lineage.



## Sequencing Information



# Initial computing setup
Before we can perform these analyses, we have to load the required software and scripts, and structure the flow of data in an organized manner. In this case, I began by organizing by dataset, as delimited by sequencing run. However, I did not detect any differences in read quality or variant detection among the separate runs, so I instead structured it based on clone category: ancestral or evolved. All of the initial analysis only requires the ancestral data, to detect existing heterozygous sites and standing variation.

All of the server to server movement is performed with rclone, while local to remote connections are through ssh. 

## Transfer fastQ files from the IGM ftp server to the Lab Google Drive


## Transfer files from Drive to TSCC



## Download and install software to TSCC
Software is installed to the home drive, as it is not temporary storage (always back up your scratch drive!).


## Upload scripts to TSCC


# Pipeline scripts and descriptions




Figure file names and descriptions
=================================

refAlelleDP-GT_underDP100
![refAlelleDP-GT_underDP100](/assets/images/refAlelleDP-GT_underDP100.pdf)
