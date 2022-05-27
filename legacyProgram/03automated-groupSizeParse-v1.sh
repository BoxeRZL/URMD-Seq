#!/bin/bash

# $laneListFile = list of lane names
# $dataLocation = folder with raw fastq files
# trim quality and flash minimum overlap are hardcoded to 30 and 5 respectivelyecho $laneListFile
trimQuality=30
flashOverlap=5

echo $laneListFile
echo $dataLocation

cat $laneListFile | while read Sample
do
	echo $Sample
	export laneName=$Sample
	python 03parseGroupSizes-v1.py $Sample $dataLocation/$Sample-$trimQuality-$flashOverlap
	cat $dataLocation/$Sample-$trimQuality-$flashOverlap/groupSizesSummary-$Sample.csv >> $dataLocation/groupSizesSummary.csv
	#cp $dataLocation/$Sample-$trimQuality-$flashOverlap/groupSizesSummary-$Sample.csv ./R3-vpcRepeat/R4-groupSizeSummaries/groupSizesSummary-$Sample.csv
done

echo "Done"
