#!/bin/bash
#PBS -j oe
#PBS -N tajimaD_test
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

module load bedtools2
module load vcf-kit 
module load bcftools


cd 2024analysis

cal_td(){
pop=$1
bcftools view -S data/population_lists/${pop}.txt data/gatk.filtered.relaxed_studywide.vcf.gz | vk tajima 10000 2000 - | sed '1d' | awk '{print $1"\t"$2+1"\t"$3"\t"$6}' > ${pop}_1based.txt
} 

export -f cal_td 
cal_td sud cha ckc her kes rhq stb ton tor vla

