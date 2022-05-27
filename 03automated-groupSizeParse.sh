#!/bin/bash

module load mafft python/3.7.0 trimmomatic

source $1
echo $laneListFile
echo $dataLocation

mkdir -p $dataLocation/$runName-groupSizeSummaries

cat $laneListFile | while read Sample
do
	echo $Sample
	export laneName=$Sample
	python 03parseGroupSizes.py $Sample $dataLocation/$Sample-$trimThreshold-$FlashThreshold
	cat $dataLocation/$Sample-$trimThreshold-$FlashThreshold/groupSizesSummary-$Sample.csv >> $dataLocation/$runName-groupSizeSummaries/$runName-groupSizes-all.csv
	cp $dataLocation/$Sample-$trimThreshold-$FlashThreshold/groupSizesSummary-$Sample.csv $dataLocation/$runName-groupSizeSummaries/groupSizesSummary-$Sample.csv
done

echo "Done"
