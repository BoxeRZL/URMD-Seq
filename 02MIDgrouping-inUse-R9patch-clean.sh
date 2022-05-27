###Dependencies: trimmomatic-0.39, flash-1.2.11, mafft/2.6
###Usage: PID-dloop-wrapper.sh LaneName /fastq/Location

module load mafft python trimmomatic

echo converting lane:  $1
echo fastq location: $2
echo trimmomatic $3 flash $4
echo move files back $6
echo saved working directory $7

#NOTE: $2 saved directory ends in '/' 
mkdir $2/$1-$3-$4/$1-{1,2,3,4,5,6,7,8,10,11,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45}
mv $2/$1-$3-$4/rawFiles/$1-$3-$4.extendedFrags.fastq $2/$1-$3-$4/

echo 'multiprocessing'
python ./02parseGroups-WTadjusted-R9patch-clean.py $2/$1-$3-$4 $1-$3-$4.extendedFrags.fastq 00cambridgeRefLong.txt m $1


homeVariable=$(pwd)
for mid in {1,2,3,4,5,6,7,8,10,11,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45}
do
	cd $homeVariable
	cd $2/$1-$3-$4/$1-$mid
	zip -mq ./alignedGroups-$1-$mid.zip ./aligned*.aln
	zip -mq ./rawGroups-$1-$mid.zip ./group*.txt
	zip -mq ./WTgroups-$1-$mid.zip ./WTgroup*.txt
	zip -mq ./WTgroupsAligned-$1-$mid.zip ./WTgroup*.aln
done
cd $SLURM_TMPDIR
echo $(pwd)
echo $1
echo $savedWd
ls
echo $6/$1-$3-$4 $savedWd/$6
cp -r $6/$1-$3-$4 $savedWd/$6
