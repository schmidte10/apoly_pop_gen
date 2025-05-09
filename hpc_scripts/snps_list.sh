#!/bin/bash
#PBS -j oe
#PBS -N snp.list
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

shopt -s expand_aliases
source /etc/profile.d/modules.sh   

cd ./2024analysis/04.OutlierAnalysis

VCF=data/gatk.filtered.relaxed_studywide.vcf.gz 
zcat $VCF | grep -v "^#" | cut -f1-4 | awk '{print $0"\t"NR}' > analysis/apoly_populations_SNPs.txt 

cd analysis/pcadapt/ 

awk '{print $2}' apoly_populations_pcaoutliers.txt > apoly_populations_pcaoutliers_numbers.txt 
awk 'FNR==NR{a[$1];next} (($4) in a)' apoly_populations_pcaoutliers_numbers.txt ../apoly_populations_SNPs.txt | cut -f4 > pcadapt_outlierSNPIDs.txt 
head pcadapt_outlierSNPIDs.txt
