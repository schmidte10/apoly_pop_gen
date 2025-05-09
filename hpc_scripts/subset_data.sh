#!/bin/bash
#PBS -j oe
#PBS -N subset_data
#PBS -l select=1:ncpus=2:mem=200gb
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

cd ./2024analysis

VCF=data/gatk.filtered.relaxed_studywide.vcf.gz
output_file=./data_subsets/"subset_${i}.vcf.gz"

for i in $(seq 0 199); do bcftools view -H $VCF | awk "NR%200==$i" | bcftools view -h $VCF | gzip -c > $output_file 

done