#!/bin/bash
#PBS -j oe
#PBS -N pstrtfltr
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

module load conda3 
shopt -s expand_aliases
source /etc/profile.d/modules.sh   

cd ./2024analysis/03.Ordinations

module load vcftools
module load bcftools
 
vcftools --gzvcf ../gatk/gatk.filtered.vcf.gz \
--max-missing 0.5 --min-meanDP 10 --max-meanDP 30 \
--minDP 3 --minGQ 40 --mac 4 --remove-filtered-geno-all --recode \
--recode-INFO-all --stdout | \
gzip > gatk.filtered.stringent_studywide.vcf.gz 

vcftools --remove-indv CTON069_S13 --gzvcf gatk.filtered.stringent_studywide.vcf.gz \
--recode --recode-INFO-all --stdout | gzip > gatk.filtered.stringent_studywide.vcf.gz


bcftools query -f '%POS\n' gatk.filtered.stringent_studywide.vcf.gz | wc -l 
bcftools query -l gatk.filtered.stringent_studywide.vcf.gz| wc -l 

vcftools --gzvcf ../gatk/gatk.filtered.vcf.gz \
--max-missing 0.5 --min-meanDP 10 --max-meanDP 30 \
--minDP 3 --minGQ 40 --mac 2 --remove-filtered-geno-all --recode \
--recode-INFO-all --stdout | \
gzip > gatk.filtered.relaxed_studywide.vcf.gz 

vcftools --remove-indv CTON069_S13 --gzvcf gatk.filtered.relaxed_studywide.vcf.gz \
--recode --recode-INFO-all --stdout | gzip > gatk.filtered.relaxed_studywide2.vcf.gz

bcftools query -f '%POS\n' gatk.filtered.relaxed_studywide.vcf.gz | wc -l 
bcftools query -l gatk.filtered.relaxed_studywide.vcf.gz| wc -l 






 

