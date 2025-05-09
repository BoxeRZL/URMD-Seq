#!/bin/bash

#SBATCH --time=0-11:59
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4000M

#SBATCH --job-name=R14-Analysis
#SBATCH --mail-user=zesli@student.ubc.ca
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

module load mafft python/3.7.0 trimmomatic
pip install biopython
#####export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
#####export MKL_NUM_THREADS=1


echo $laneName
echo $WTadjustedFlag

# Guaranteed to use Python 3.X
python3 04writeMutations-multiprocessing-R14patch-clean-debug.py $dataLocation/$laneName-$trimThreshold-$FlashThreshold 5 0.75 00cambridgeRefLong.txt $dataLocation $outputFileName $laneName- $WTadjustedFlag


echo "Done"
