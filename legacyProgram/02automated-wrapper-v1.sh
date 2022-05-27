#!/bin/bash

# $laneListFile = list of lane names
# $dataLocation = folder with raw fastq files
# trim quality and flash minimum overlap are hardcoded to 30 and 5 respectivelyecho $laneListFile
trimQuality=30
flashOverlap=5

cat $laneListFile | while read Sample
do
	echo $Sample
	export laneName=$Sample
	bash 02MIDgrouping-v1.sh $laneName $dataLocation $trimQuality $flashOverlap -m

done

echo "Done"
