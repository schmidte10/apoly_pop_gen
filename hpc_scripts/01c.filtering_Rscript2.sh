#!/bin/bash
#PBS -j oe
#PBS -N filter_snps_plotR2
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

module load conda3 
shopt -s expand_aliases
source /etc/profile.d/modules.sh   

cd ./2024analysis/gatk

module load R 

R -f filtering_Rscript2.R
