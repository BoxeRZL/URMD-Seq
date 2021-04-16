# dloop-PID: Mitochondrial d-loop Primer ID Sequencing
--------------

This repository is a self contained analysis pipeline for targetted mitochondrial d-loop sequencing using a unique primer barcoding technique.  This pipeline takes Illumina 2X300 paired reads fastq files, spanning the d-loop region of the mitochondrial DNA (M:16559-279).  

Requirements: Python/3.7.0 <br>
Dependencies: Biopython/1.75, Mafft/2.6, Trimmomatic/0.39, FLASH/1.2.11

Data set up requirements:
A list of lane names should be in a text file with one name per line, these refer to the identifier preceeding `\_R\*\_001.fastq` in the raw fastq files. This allows for automated running of all lanes.

----------------------

1. Initial preprocessing trims Illumina adapters and filters based on both sliding window and average quality scores, combines paired reads with a minimum overlap. The percentage of reads surviving trimming and read combining can be found in the data folder, labeled as 'preprocessingQC.csv'

  `export dataLocation=./R4-adam
  export laneListFile=R4-laneList.txt`

  `bash 01preprocessing-QC-automated.sh -export=dataLocation -export=laneListFile`

  `dependencies
01PID-dloop-preprocessing-QC.sh
01preprocessing-QC-automated.sh`

<br><br>

2. Barcoded reads are grouped together and aligned for subsequent analysis step, all barcodes with only one read present are not grouped and are appended to the file 'group0.txt'.  All sequences without a recognizable ID can be found in the file 'mutated.txt'

 `export dataLocation=./R4-adam
  export laneListFile=R4-laneList-1.txt`

 `bash 02automated-wrapper-v1.sh -export=dataLocation -export=laneListFile`

  `dependencies
02MIDgrouping-v1.sh
02parseGroups-v1.py`

By default python multiprocessing is turned on (using 'multiprocessing' library)
<br><br>

3. The number of constituent reads for each strand consensus are summarized in the file 'groupSizesSummary.csv'. This is used to ensure there was sufficient amplification and that there is an even distribution of reads across strand consensuses.

 `export dataLocation=./R4-adam
  export laneListFile=R4-laneList-1.txt`

 `bash 03automated-groupSizeParse-v1.sh -export=dataLocation -export=laneListFile`

  `dependencies
03parseGroupSizes-v1.py`

<br><br>

4. Mutation calling parses low level mutations, only accepting valid sequences if the number of reads in a barcoded group is above a threshold (default=5), and if the consensus within the group is above a threshold (default=75%).  The cambridge reference sequence (NC_012920) can be obtained [here](https://www.ncbi.nlm.nih.gov/nuccore/251831106), for this purpose positions 16559-279 was used

 `export dataLocation=./R4-adam
export laneListFile=R4-laneList-1.txt
export outputFileName=Example-outputFile-`

Default reversion of ambiguous bases to wild type allele
Default minimum constituent read group size and mutation frequency threshold are set to 5 and 0.75 respectively

`bash 04automated-writeMutations-v1.sh -export=laneListFile -export=dataLocation -export=outputFileName`

`dependencies
04writeMutations-v1.py`




--------------------

Analysis Schematic

![](https://github.com/racheldunn/dloop-PID/blob/master/extraFiles/dloop-PID-schematic.png "dloop-PID-schematic")

\*step 1 outlined in orange, step 2 outlined in blue\*

<br><br>
--------------------
To do: 
1. update contamination filter
2. update slurm enabled version
