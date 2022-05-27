#!/bin/bash

module load mafft python/3.7.0 trimmomatic
source $1
echo $laneListFile

cat $laneListFile | while read Sample
do
	echo $Sample
	export trimThreshold=$trimThreshold
	export FlashThreshold=$FlashThreshold
	export laneName=$Sample
	export savedWd=$(pwd)
	export dataLocation=$dataLocation

	sbatch 02sbatch-MIDgrouping-R9patch-clean.sh -export=laneName -export=dataLocation -export=savedWd -export=trimThreshold -export=FlashThreshold
done

echo "Done"
