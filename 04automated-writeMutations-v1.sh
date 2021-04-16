#!/bin/bash

# $laneListFile = list of lane names
# $dataLocation = folder with raw fastq files
# $outputFileName = processed output will be printed to this file (the date and file type are appended to what is specified here)
# trim quality and flash minimum overlap are hardcoded to 30 and 5 respectively
trimQuality=30
flashOverlap=5

# minimum constituent read group size and mutation frequency threshold are hardcoded to 5 and 0.75 respectively
groupSize=5
frequencyThreshold=0.75

#WTadjustFlag default is on, see details in manuscript
WTadjustedFlag=True

cat $laneListFile | while read Sample
do
	echo $Sample
	export laneName=$Sample
	#sbatch 04sbatch-writeMutations.sh -export=laneName -export=dataLocation -export=outputFileName -export=WTadjustedFlag
	python 04writeMutations-v1.py $dataLocation/$laneName-$trimQuality-$flashOverlap $groupSize $frequencyThreshold 00cambridgeRefLong.txt $dataLocation $outputFileName $laneName- $WTadjustedFlag

done

echo "Done"
