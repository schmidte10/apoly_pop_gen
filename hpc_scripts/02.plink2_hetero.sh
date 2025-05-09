#!/bin/bash
#PBS -j oe
#PBS -N heterozygosity
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

module load plink2 

cd ./2024analysis
plink2 --vcf gatk/gatk.filteredv3.vcf.gz --allow-extra-chr --het --bad-freqs