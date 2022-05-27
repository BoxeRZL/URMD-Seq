# dloop-PID: Mitochondrial d-loop Primer ID Sequencing
--------------

This repository is a self contained analysis pipeline for targetted mitochondrial d-loop sequencing using a unique primer barcoding technique.  This pipeline takes Illumina 2X300 paired reads fastq files, spanning the d-loop region of the mitochondrial DNA (M:16559-279).  

<br>
Dependencies: Python/3.7.0 Biopython/1.75, Mafft/2.6, Trimmomatic/0.39, FLASH/1.2.11
<br><br>
Data set up requirements: <br>
A list of lane names should be in a text file with one name per line, these refer to the identifier preceeding `\_R\*\_001.fastq` in the raw fastq files. This allows for automated analysis of all lanes.
<br><br>
All analyses are performed using Compute Canada with modules:
`StdEnv/2018.3`
`mafft python/3.7 trimmomatic`
<br><br>

----------------------


Configuration file:

`runName=RX`\*<br>
`dataLocation=RX_sequencingData`<br>
`laneListFile=RX-laneList.txt`\*\*<br>
`trimThreshold=30`\*\*\*<br>
`FlashThreshold=5`****<br>
`outputFileName=RX-mutations`

\*Where ‘X’ is the run number (ie. R5), and dataLocation is the folder where the fastq files are (ie R5_sequencingData)<br>
\*\*text file with lane names<br>
\*\*\*Trim threshold is the minimum quality score that is acceptable during trimming, 30 is the default, may be increased for more stringency or decreased if the run data is of low quality (Trimmomatic, http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf )<br>
\*\*\*\*Flash threshold is the minimum overlap required between the forward and reverse reads, for R5 and R6 it was 5 and for R4 it was 30 (FLASH, http://ccb.jhu.edu/software/FLASH/MANUAL )<br>

----------------------


1. Initial preprocessing trims Illumina adapters and filters based on both sliding window and average quality scores, combines paired reads with a minimum overlap. The percentage of reads surviving trimming and read combining can be found in the data folder, labeled as 'preprocessingQC.csv'


  `bash 01preprocessing-QC-automated.sh ./config.fileRX`
<br><br>
dependencies<br>
`01PID-dloop-preprocessing-QC.sh`<br>
`01preprocessing-QC-automated.sh`<br>

<br>

2. Barcoded reads are grouped together and aligned for subsequent analysis step, all barcodes with only one read present are not grouped and are appended to the file 'group0.txt'.  All sequences without a recognizable ID can be found in the file 'mutated.txt'
sbatch parameters specified in `02sbatch-MIDgrouping-R9patch-clean.sh` edit as necessary

  `bash 02automated-wrapper-R9patch-clean.sh ./config.fileRX`
<br><br>
dependencies<br>
`02sbatch-MIDgrouping-R9patch-clean.sh`<br>
`02MIDgrouping-inUse-R9patch-clean.py`<br>
`02parseGroups-WTadjusted-R9patch-clean.py`<br><br>

By default python multiprocessing is turned on (using 'multiprocessing' library)
<br><br>

3. The number of constituent reads for each strand consensus are summarized in the file 'groupSizesSummary.csv'. This is used to ensure there was sufficient amplification and that there is an even distribution of reads across strand consensuses.

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


