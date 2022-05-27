#!/bin/bash

###Dependencies: trimmomatic-0.39, flash-1.2.11, mafft/2.6
###Usage: PID-dloop-wrapper.sh LaneName /fastq/Location

module load mafft python/3.7.0 trimmomatic

echo converting lane:  $1
echo fastq location: $2
echo trimmomatic $3 flash $4

mkdir $2/$1-$3-$4
mkdir $2/$1-$3-$4/rawFiles
#cp $2/$1_R1_001.fastq $2/$1-$3-$4/rawFiles/

#RP=$(awk '{if(NR==4) print $0}' ./R2-Jan2020/R2-30-5-preprocessingQC.txt | awk '{print $4}')
#touch $2/$1-$3-$4-preprocessingQC-R1-updatedScript.csv

java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.36.jar PE -phred33 $2/$1_R1_001.fastq $2/$1_R2_001.fastq -baseout $2/$1-$3-$4/rawFiles/$1-$3-$4_trimmed.fastq SLIDINGWINDOW:4:$3 AVGQUAL:$3 |& tee $2/$1-$3-$4-$6
~/projects/def-hecote/dloop_shared/FLASH-1.2.11-Linux-x86_64/flash -d $2/$1-$3-$4/rawFiles -o $1-$3-$4 $2/$1-$3-$4/rawFiles/$1-$3-$4_trimmed_1P.fastq $2/$1-$3-$4/rawFiles/$1-$3-$4_trimmed_2P.fastq -m $4 -M 360 |& tee -a $2/$1-$3-$4-$6

#java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.36.jar PE -phred33 $2/$1/rawFiles/$1_R1_001.fastq $2/$1/rawFiles/$1_R2_001.fastq -baseout $2/$1/rawFiles/$1_trimmed.fastq SLIDINGWINDOW:4:30 AVGQUAL:30
#/home/radunn/projects/def-hecote/radunn/myPackages/FLASH-1.2.11-Linux-x86_64/flash -m 100 -M 400 -d $2/$1/rawFiles -o $1 $2/$1/rawFiles/$1_trimmed_1P.fastq $2/$1/rawFiles/$1_trimmed_2P.fastq



#mkdir $2/$1-$3-$4/$1-{1,3,4,5,6,7,10,13,14,15,16,17,18}
#mv $2/$1-$3-$4/rawFiles/$1-$3-$4.extendedFrags.fastq $2/$1-$3-$4/

#if [ $# -eq 5 ]; then
#	echo 'multiprocessing'
#	python ./grpParse-WTadjusted-dec2020.py $2/$1-$3-$4 $1-$3-$4.extendedFrags.fastq cambridgeRefLong.txt m $1
#else
#	echo 'singleprocessing'
#	python ./grpParse-WTadjusted-dec2020.py $2/$1-$3-$4 $1-$3-$4.extendedFrags.fastq cambridgeRefLong.txt s $1
#fi

#cd $2/$1-$3-$4/$1-1
#zip -mq ./alignedGroups-$1-1.zip ./aligned*.aln
#zip -mq ./rawGroups-$1-1.zip ./group*.txt
#zip -mq ./WTgroups-$1-1.zip ./WTgroup*.aln

#for mid in {3,4,5,6,7,10,13,14,15,16,17,18}
#do
#	echo $1-$mid
#	cd ../$1-$mid
#	zip -mq ./alignedGroups-$1-$mid.zip ./aligned*.aln
#	zip -mq ./rawGroups-$1-$mid.zip ./group*.txt
#	zip -mq ./WTgroups-$1-$mid.zip ./WTgroup*.txt
#done

#cd ../1000Apg-1_S3/1
#zip -mq ./alignedGroups-$1-$mid.zip ./aligned*.aln
#zip -mq ./rawGroups-$1-$mid.zip ./group*.txt

#for mid in {4,5,6,7,13,15,16,17}
#do
#	echo $mid
#	cd ../$mid
#	pwd
#	zip -mq ./alignedGroups$mid.zip ./aligned*.aln
#	zip -mq ./rawGroups$mid.zip ./group*.txt
#done
