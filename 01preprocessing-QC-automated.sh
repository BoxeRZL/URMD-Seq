#!/bin/bash

inputSuffix='preprocessingQC.txt'

module load trimmomatic mafft python/3.7.0

source $1
echo $trimThreshold
echo $FlashThreshold

echo 'SampleName ReadPairs BothSurviving BothSurvivingPercent OnlyForwardSurviving OnlyForwardPercent OnlyReverseSurviving OnlyReversePercent DiscardedReads DiscardedPercent flashTotalPairs flashCombinedPairs flashUncombinedPairs flashPercentCombined' |& tee $dataLocation/$runName-$trimThreshold-$FlashThreshold-$inputSuffix

cat $laneListFile | while read Sample
do
	bash 01PID-dloop-preprocessing-QC.sh $Sample $dataLocation $trimThreshold $FlashThreshold -m $inputSuffix

	RP=$(grep -n 'Input Read Pairs:' $dataLocation/$Sample-$trimThreshold-$FlashThreshold-$inputSuffix | awk '{print $4}')
	BSR=$(grep -n 'Input Read Pairs:' $dataLocation/$Sample-$trimThreshold-$FlashThreshold-$inputSuffix | awk '{print $7}')
	BSP=$(grep -n 'Input Read Pairs:' $dataLocation/$Sample-$trimThreshold-$FlashThreshold-$inputSuffix | awk '{print $8}' | awk '{print substr($0, 2, length($0) - 2)}')

	FOR=$(grep -n 'Input Read Pairs:' $dataLocation/$Sample-$trimThreshold-$FlashThreshold-$inputSuffix | awk '{print $12}')
	FOS=$(grep -n 'Input Read Pairs:' $dataLocation/$Sample-$trimThreshold-$FlashThreshold-$inputSuffix | awk '{print $13}' | awk '{print substr($0, 2, length($0) - 2)}')

	ROR=$(grep -n 'Input Read Pairs:' $dataLocation/$Sample-$trimThreshold-$FlashThreshold-$inputSuffix | awk '{print $17}')
	ROS=$(grep -n 'Input Read Pairs:' $dataLocation/$Sample-$trimThreshold-$FlashThreshold-$inputSuffix | awk '{print $18}' | awk '{print substr($0, 2, length($0) - 2)}')

	DR=$(grep -n 'Input Read Pairs:' $dataLocation/$Sample-$trimThreshold-$FlashThreshold-$inputSuffix | awk '{print $20}')
	DS=$(grep -n 'Input Read Pairs:' $dataLocation/$Sample-$trimThreshold-$FlashThreshold-$inputSuffix | awk '{print $21}' | awk '{print substr($0, 2, length($0) - 2)}')

	TPR=$(grep -n 'Total pairs:' $dataLocation/$Sample-$trimThreshold-$FlashThreshold-$inputSuffix | awk '{print $4}')
	CPR=$(grep -n 'Combined pairs' $dataLocation/$Sample-$trimThreshold-$FlashThreshold-$inputSuffix | awk '{print $4}')
	UPR=$(grep -n 'Uncombined pairs' $dataLocation/$Sample-$trimThreshold-$FlashThreshold-$inputSuffix | awk '{print $4}')
	CS=$(grep -n 'Percent combined' $dataLocation/$Sample-$trimThreshold-$FlashThreshold-$inputSuffix | awk '{print $4}')
	echo $Sample $RP $BSR $BSP $FOR $FOS $ROR $ROS $DR $DS $TPR $CPR $UPR $CS |& tee -a $dataLocation/$runName-$trimThreshold-$FlashThreshold-$inputSuffix
	rm $dataLocation/$Sample-$trimThreshold-$FlashThreshold-$inputSuffix
done
