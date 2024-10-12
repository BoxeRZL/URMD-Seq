#!/bin/bash

module load mafft python/3.7.0
source $1

echo $laneListFile
echo $WTadjustFlag

cat $laneListFile | while read Sample
do
	echo $Sample
	export laneName=$Sample
	export WTadjustedFlag=True
	export dataLocation=$dataLocation
	export outputFileName=$outputFileName
	export trimThreshold=$trimThreshold
	export FlashThreshold=$FlashThreshold
	sbatch 04sbatch-writeMutations-R14patch-clean.sh -export=laneName -export=dataLocation -export=outputFileName -export=WTadjustedFlag -export=trimThreshold -export=FlashThreshold
done

echo "Done"
