#!/bin/bash
#PBS -j oe
#PBS -N PCAdapt
#PBS -l select=1:ncpus=2:mem=100gb
#PBS -l walltime=05:00:00
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

cd 2024analysis

module load vcftools 
module load plink

# VCF to PLINK
vcftools --vcf data/gatk.filtered.relaxed_studywide.vcf.gz --plink --out data/gatk.filtered.relaxed_studywide.plink 

# convert PLINK to BED
plink --file data/gatk.filtered.relaxed_studywide.plink --make-bed --noweb --out apoly_populations

module load R 

R -f Rscripts/PCAdapt_script1.R

