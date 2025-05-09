#!/bin/bash
#PBS -j oe
#PBS -N Bayescan
#PBS -l select=1:ncpus=2:mem=100gb
#PBS -l walltime=30:00:00
#PBS -m ae
#PBS -M elliott.schmidt@my.jcu.edu.au

echo "------------------------------------------------------"
echo "PBS: Submitted to $PBS_QUEUE@$PBS_O_HOST"
echo "PBS: Working directory is $PBS_O_WORKDIR"
echo "PBS: Job identifier is $PBS_JOBID"
echo "PBS: Job name is $PBS_JOBNAME$no"
echo "------------------------------------------------------"

shopt -s expand_aliases
source /etc/profile.d/modules.sh   

cd ./2024analysis/04.OutlierAnalysis

module load bayescan 

bayescan ./analysis/bayescan/apoly_populations2.bs -od ./analysis/bayescan -threads 2 -n 5000  -thin 10 -nbp 20 -pilot 5000 -burn 50000 -pr_odds 10
