#!/bin/bash
#PBS -j oe
#PBS -N plink2pca
#PBS -l select=1:ncpus=2:mem=120gb
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

cd ./2024analysis

module load plink2

plink2 --vcf ./gatk/gatk.filteredv3.vcf.gz --make-pgen --allow-extra-chr --out ./03.Ordinations/gatk.filteredv3
plink2 --pfile ./03.Ordinations/gatk.filteredv3 --allow-extra-chr --pca --out ./03.Ordinations/gatk.filteredv3.pca





 

