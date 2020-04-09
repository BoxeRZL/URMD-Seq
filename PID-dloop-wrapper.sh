#!/bin/bash

###Dependencies: trimmomatic-0.39, flash-1.2.11, mafft/2.6
###Usage: PID-dloop-wrapper.sh LaneName /fastq/Location

echo converting lane:  $1
echo fastq location: $2

mkdir $2/$1
mkdir $2/$1/rawFiles
mv $2/$1_R* $2/$1/rawFiles/

java -jar /Users/racheldunn/masters/packagesLoaded/Trimmomatic-0.39/trimmomatic-0.39.jar PE -phred33 $2/$1/rawFiles/$1_R1_001.fastq $2/$1/rawFiles/$1_R2_001.fastq -baseout $2/$1/rawFiles/$1_trimmed.fastq SLIDINGWINDOW:4:30 AVGQUAL:30
/Users/racheldunn/masters/packagesLoaded/FLASH-1.2.11/flash -m 100 -M 400 -d $2/$1/rawFiles -o $1 $2/$1/rawFiles/$1_trimmed_1P.fastq $2/$1/rawFiles/$1_trimmed_2P.fastq 

mkdir $2/$1/{1,4,5,6,7,13,15,16,17}
mv $2/$1/rawFiles/$1.extendedFrags.fastq $2/$1/

if [ $# -eq 3 ]; then
	echo 'multiprocessing'
	python ./grpParseClean-2020.py $2/$1 $1.extendedFrags.fastq cambridgeRefLong.txt m
else
	echo 'singleprocessing'
	python ./grpParseClean-2020.py $2/$1 $1.extendedFrags.fastq cambridgeRefLong.txt s
fi

#cd $2/$1/1
#zip -mq ./alignedGroups-$1-$mid.zip ./aligned*.aln
#zip -mq ./rawGroups-$1-$mid.zip ./group*.txt

for mid in {1,4,5,6,7,13,15,16,17}
do
	echo $mid
	cd $2/$1/$mid
	zip -mq ./alignedGroups-$1-$mid.zip ./aligned*.aln
	zip -mq ./rawGroups-$1-$mid.zip ./group*.txt
done

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
