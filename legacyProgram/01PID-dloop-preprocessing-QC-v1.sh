#!/bin/bash

###Dependencies: trimmomatic-0.39, flash-1.2.11, mafft/2.6
# bash 01PID-dloop-preprocessing-QC-v1.sh $Sample $dataLocation $trimQuality $flashOverlap -m

echo converting lane:  $1
echo fastq location: $2
echo trimmomatic $3 flash $4

mkdir $2/$1-$3-$4
mkdir $2/$1-$3-$4/rawFiles

java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.36.jar PE -phred33 $2/$1_R1_001.fastq $2/$1_R2_001.fastq -baseout $2/$1-$3-$4/rawFiles/$1-$3-$4_trimmed.fastq SLIDINGWINDOW:4:$3 AVGQUAL:$3 |& tee $2/$1-$3-$4-preprocessingQC.txt
~/projects/def-hecote/dloop_shared/FLASH-1.2.11-Linux-x86_64/flash -d $2/$1-$3-$4/rawFiles -o $1-$3-$4 $2/$1-$3-$4/rawFiles/$1-$3-$4_trimmed_1P.fastq $2/$1-$3-$4/rawFiles/$1-$3-$4_trimmed_2P.fastq -m $4 -M 400 |& tee -a $2/$1-$3-$4-preprocessingQC.txt
