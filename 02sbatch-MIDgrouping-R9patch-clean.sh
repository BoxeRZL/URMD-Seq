#!/bin/bash

#SBATCH --time=0-71:59
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=18000M

#SBATCH --job-name=R7Analysis
#SBATCH --mail-user=racheld3141@gmail.com
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

module load mafft python/3.7.0 trimmomatic
pip install biopython
#####export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
#####export MKL_NUM_THREADS=1

echo $laneListFile

cp -r $dataLocation $SLURM_TMPDIR
bash 02MIDgrouping-inUse-R9patch-clean.sh $laneName $SLURM_TMPDIR/$dataLocation 30 5 -m $dataLocation $savedWd

echo "Done"
