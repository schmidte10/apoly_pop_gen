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

module load bcftools 
module load tabix
module load plink

bcftools view -m2 -M2 -v snps data/gatk.filtered.relaxed_studywide.bicolor.vcf.gz -Oz -o data/gatk.filtered.relaxed_studywide.bicolor.biallelic.vcf.gz 

tabix -p vcf data/gatk.filtered.relaxed_studywide.bicolor.biallelic.vcf.gz

plink2 --vcf data/gatk.filtered.relaxed_studywide.bicolor.biallelic.vcf.gz --make-bed --allow-extra-chr --out data/gatk.filtered.relaxed_studywide.bicolor.biallelic.test.vcf.gz
module load R 

R -f Rscripts/PCAdapt_script1.R

