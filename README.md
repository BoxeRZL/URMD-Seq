# UMI-DloopSeq: Mitochondrial d-loop UMI based Sequencing

--------------

Contact: Zeshuo Li (lizeshuo99@126.com)

--------------

This repository is a self contained analysis pipeline for targetted mitochondrial d-loop sequencing using a unique primer barcoding technique.  This pipeline takes Illumina 2X300 paired reads fastq files, spanning the d-loop region of the mitochondrial DNA (M:16559-279).  

<br>
Dependencies: Python/3.7.0, Biopython/1.75, Mafft/2.6, Trimmomatic/0.39, FLASH/1.2.11
<br><br>

Data formatting requirements: <br>

All data should be in raw fastq format with the suffix `\_RX\_001.fastq` , for automated analysis a text file is created containing a list of all lane names (preceeding `\_RX\_001.fastq`)
<br><br>
All analyses are performed using Digital Research Alliance Canada (used called Compute Canada), with slurm and modules:
`StdEnv/2018.3`
`mafft python/3.7 trimmomatic`
<br><br>

----------------------


Configuration file:

`runName=RX`\*<br>
`dataLocation=RX_sequencingData`\*\*<br>
`laneListFile=RX-laneList.txt`\*\*\*<br>
`trimThreshold=30`\*\*\*\*<br>
`FlashThreshold=5`*****<br>
`outputFileName=RX-mutations`

\*‘RX’ is the name of the run (i.e. R11; edit it to the name you need)<br>
\*\*Data location is the folder where row FASTQ files are (i.e. R11-sequencingData)<br>
\*\*\*Lane list file is a test files containing a list of all lane names (I.e. R11-laneList.txt). It is for automated analysis to proceed all FASTQ files (i.e. “_R1_001.fastq” and “_R2_001.fastq”)<br>
\*\*\*\*Trim threshold is the minimum quality score that is acceptable during trimming. 30 is the default. The threshold may be increased for more stringency or decreased if the sequencing data is of low quality (see the documentation of Trimmomatic for details, http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf )<br>
\*\*\*\*\*Flash threshold is the minimum length of overlap required between the forward and reverse reads. 10 is the default. Various flash thresholds were tested by using our data and 5 is the optimized threshold. It may be edited if necessary.For R5 and R6 it was 5 and for R4 it was 30 (FLASH, http://ccb.jhu.edu/software/FLASH/MANUAL )<br>

----------------------


1.Trim sequencing adapters of Illumina and filter out reads with poor-quality bases and then combine paired reads on a minimum length of overlap. The default settings are specified in the configuration file ‘config.fileRX’. The reads surviving of trimming and combining are counts in percentages for quality control and can be found in an output file labeled as ‘preprocessingQC.csv’ in the data location folder (see Box2 for interpretations of quality control parameters).


  `bash 01preprocessing-QC-automated.sh ./config.fileRX`
<br><br>
dependencies<br>
`01PID-dloop-preprocessing-QC.sh`<br>
`01preprocessing-QC-automated.sh`<br>

<br>

2.Group barcoded reads with same UMI into a family. All UMI barcodes with the presence of less than two reads are not grouped and are appended to an output file labeled as ‘ungrouped.txt’. All reads without a recognizable UMI are listed in an output file labeled as ‘mutated.txt’.  An initial alignment is processed between reads within a family to generate a family consensus which is appended to an output file labeled as ‘Consensus.aln’. The number of reads in each family is calculated into an output file labeled as 'groupsize.csv' for subsequent analysis steps. <br>
sbatch parameters specified in `02sbatch-MIDgrouping-R9patch-clean.sh` edit as necessary

  `bash 02automated-wrapper-R9patch-clean.sh ./config.fileRX`
<br><br>
dependencies<br>
`02sbatch-MIDgrouping-R9patch-clean.sh`<br>
`02MIDgrouping-inUse-R9patch-clean.py`<br>
`02parseGroups-WTadjusted-R9patch-clean.py`<br><br>

By default python multiprocessing is turned on (using 'multiprocessing' library)
<br><br>

3.Summarize the number of constituent reads among all families in an output file labeled as ‘groupSizesSummary.csv’. This is used to indicate whether there was sufficient amplification for quality control.

  `bash 03automated-groupSizeParse.sh ./config.fileRX`
<br><br>
dependencies<br>
`03parseGroupSizes.py`


<br><br>

4. Mutation calling parses low level mutations, only accepting valid sequences if the number of reads in a barcoded group is above a threshold (default=5), and if the consensus within the group is above a threshold (default=75%).  The cambridge reference sequence (NC_012920) can be obtained [here](https://www.ncbi.nlm.nih.gov/nuccore/251831106), for this purpose positions 16559-279 was used

  `bash 04automated-writeMutations-clean.sh ./config.fileRX`<br><br>

dependencies<br>
`04sbatch-writeMutations-clean.sh`
`04writeMutations-multiprocessing-clean.py`
<br><br>

--------------------


