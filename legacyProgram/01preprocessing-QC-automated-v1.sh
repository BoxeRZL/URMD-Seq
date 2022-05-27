#!/bin/bash

# $laneListFile = list of lane names
# $dataLocation = folder with raw fastq files
# trim quality and flash minimum overlap are hardcoded to 30 and 5 respectively
trimQuality=30
flashOverlap=5

inputSufix=$trimQuality-$flashOverlap-preprocessingQC.csv

echo 'SampleName ReadPairs BothSurviving BothSurvivingPercent OnlyForwardSurviving OnlyForwardPercent OnlyReverseSurviving OnlyReversePercent DiscardedReads DiscardedPercent flashTotalPairs flashCombinedPairs flashUncombinedPairs flashPercentCombined' |& tee $dataLocation/$inputSufix

cat $laneListFile | while read Sample
do
	bash 01PID-dloop-preprocessing-QC-v1.sh $Sample $dataLocation $trimQuality $flashOverlap -m

	RP=$(awk '{if(NR==4) print $0}' $dataLocation/$Sample-$inputSufix | awk '{print $4}')
	BSR=$(awk '{if(NR==4) print $0}' $dataLocation/$Sample-$inputSufix | awk '{print $7}')
	BSP=$(awk '{if(NR==4) print $0}' $dataLocation/$Sample-$inputSufix | awk '{print $8}' | awk '{print substr($0, 2, length($0) - 2)}')

	FOR=$(awk '{if(NR==4) print $0}' $dataLocation/$Sample-$inputSufix | awk '{print $12}')
	FOS=$(awk '{if(NR==4) print $0}' $dataLocation/$Sample-$inputSufix | awk '{print $13}' | awk '{print substr($0, 2, length($0) - 2)}')

	ROR=$(awk '{if(NR==4) print $0}' $dataLocation/$Sample-$inputSufix | awk '{print $17}')
	ROS=$(awk '{if(NR==4) print $0}' $dataLocation/$Sample-$inputSufix | awk '{print $18}' | awk '{print substr($0, 2, length($0) - 2)}')

	DR=$(awk '{if(NR==4) print $0}' $dataLocation/$Sample-$inputSufix | awk '{print $20}')
	DS=$(awk '{if(NR==4) print $0}' $dataLocation/$Sample-$inputSufix | awk '{print $21}' | awk '{print substr($0, 2, length($0) - 2)}')

	TPR=$(grep -n 'Total pairs:' $dataLocation/$Sample-$inputSufix | awk '{print $4}')
	CPR=$(grep -n 'Combined pairs' $dataLocation/$Sample-$inputSufix | awk '{print $4}')
	UPR=$(grep -n 'Uncombined pairs' $dataLocation/$Sample-$inputSufix | awk '{print $4}')
	CS=$(grep -n 'Percent combined' $dataLocation/$Sample-$inputSufix | awk '{print $4}')
	echo $Sample $RP $BSR $BSP $FOR $FOS $ROR $ROS $DR $DS $TPR $CPR $UPR $CS |& tee -a $dataLocation/$inputSufix
	rm $dataLocation/$Sample$inputSufix
done
