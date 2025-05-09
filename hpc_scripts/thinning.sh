#!/bin/bash
#PBS -j oe
#PBS -N thinning
#PBS -l select=1:ncpus=2:mem=100gb
#PBS -l walltime=10:00:00
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

module load bcftools 

cd 2024analysis/data

bcftools view gatk.filtered.relaxed_studywide.vcf.gz | awk 'NR % 200 == 0 || /^#/' | bcftools view -o thinned_output.vcf.gz