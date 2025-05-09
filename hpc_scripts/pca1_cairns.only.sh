#!/bin/bash
#PBS -j oe
#PBS -N pca_cairns
#PBS -l select=1:ncpus=1:mem=100gb
#PBS -l walltime=30:00:00
#PBS -m ae
#PBS -M elliott.schmidt@my.jcu.edu.au

echo "------------------------------------------------------"
echo "PBS: Submitted to $PBS_QUEUE@$PBS_O_HOST"
echo "PBS: Working directory is $PBS_O_WORKDIR"
echo "PBS: Job identifier is $PBS_JOBID"
echo "PBS: Job name is $PBS_JOBNAME$no"
echo "------------------------------------------------------"

module load anaconda3
source $CONDA_PROF/conda.sh
conda activate R-4.0.3

R -f pca_cairns.only.R
