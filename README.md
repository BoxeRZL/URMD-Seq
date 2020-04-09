# dloop-PID: Mitochondrial d-loop Primer ID Sequencing
--------------

This repository is a self contained analysis pipeline for targetted mitochondrial d-loop sequencing using a unique primer barcoding technique.  This pipeline takes Illumina 2X300 paired reads fastq files, spanning the d-loop region of the mitochondrial DNA (M:16559-279).  

Requirements: Python/3.7 <br>
Dependencies: Biopython/1.75, Mafft/2.6, Trimmomatic/0.39, FLASH/1.2.11

----------------------

1. Initial preprocessing trims illumina adapters and filters based on both sliding window and average quality scores, combines paired reads with a minimum overlap.  Barcoded reads are grouped together and aligned for subsequent analysis step, all barcodes with only one read present are not grouped and are considered 'group0'.  All sequences without a recognizable ID can be found in the file 'mutated.txt'

  `Usage: PID-dloop-wrapper.sh LaneName /fastq/Location -m`

LaneName refers to the identifier preceeding `\_R\*\_001.fastq` in the raw fastq files<br>
-m turns on multiprocessing, default is off 

<br><br>

2. Mutation calling parses low level mutations, only accepting valid sequences if the number of reads in a barcoded group is above a threshold (default=5), and if the consensus within the group is above a threshold (default=75%).  The cambridge reference sequence (NC_012920) can be obtained [here](https://www.ncbi.nlm.nih.gov/nuccore/251831106), for this purpose positions 16559-279 was used

  `Usage: python ./writeMutations-2020.py /sorted/mid/folders/root groupSize(5) consensus(0.75) cambridgeRefDloop.txt`




--------------------

Analysis Schematic

![](https://github.com/racheldunn/dloop-PID/blob/master/extraFiles/dloop-PID-schematic.png "dloop-PID-schematic")

\*step 1 outlined in orange, step 2 outlined in blue\*
