# Dependencies: trimmomatic-0.39, flash-1.2.11, mafft/2.6, Biopython/1.75
# 02MIDgrouping-v1.sh $laneName $dataLocation $trimQuality $flashOverlap -m

echo converting lane:  $1
echo fastq location: $2
echo trimmomatic $3 flash $4

mkdir $2/$1-$3-$4/$1-{1,3,4,5,6,7,10,13,14,15,16,17,18}
mv $2/$1-$3-$4/rawFiles/$1-$3-$4.extendedFrags.fastq $2/$1-$3-$4/

if [ $# -eq 5 ]; then
	echo 'multiprocessing'
	python ./02parseGroups-v1.py $2/$1-$3-$4 $1-$3-$4.extendedFrags.fastq 00cambridgeRefLong.txt m $1
else
	echo 'singleprocessing'
	python ./02parseGroups-v1.py $2/$1-$3-$4 $1-$3-$4.extendedFrags.fastq 00cambridgeRefLong.txt s $1
fi

cd $2/$1-$3-$4/$1-1
zip -mq ./alignedGroups-$1-1.zip ./aligned*.aln
zip -mq ./rawGroups-$1-1.zip ./group*.txt
zip -mq ./WTgroups-$1-1.zip ./WTgroup*.txt
zip -mq ./WTgroupsAligned-$1-1.zip ./WTgroup*.aln

for mid in {3,4,5,6,7,10,13,14,15,16,17,18}
do
	echo $1-$mid
	cd ../$1-$mid
	zip -mq ./alignedGroups-$1-$mid.zip ./aligned*.aln
	zip -mq ./rawGroups-$1-$mid.zip ./group*.txt
	zip -mq ./WTgroups-$1-$mid.zip ./WTgroup*.txt
	zip -mq ./WTgroupsAligned-$1-$mid.zip ./WTgroup*.aln
done
