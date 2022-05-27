#!/bin/bash

#SBATCH --time=0-11:59
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4000M

#SBATCH --job-name=R5-Jan2021Analysis
#SBATCH --mail-user=racheld3141@gmail.com
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

module load mafft python/3.7.0 trimmomatic
pip install biopython
#####export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
#####export MKL_NUM_THREADS=1

echo $laneName
echo $WTadjustedFlag

python 04writeMutations-multiprocessing-clean.py $dataLocation/$laneName-$trimThreshold-$FlashThreshold 5 0.75 00cambridgeRefLong.txt $dataLocation $outputFileName $laneName- $WTadjustedFlag

echo "Done"
