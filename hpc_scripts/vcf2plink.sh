#!/bin/bash
#PBS -j oe
#PBS -N vcf2plink
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

cd ./2024analysis/04.OutlierAnalysis

module load vcftools
module load plink

VCF_TEST=data/freebayes_f3_chr1 
VCF=data/gatk.filtered.relaxed_studywide.vcf.gz 

vcftools --gzvcf $VCF --plink --out data/apoly_populations.plink 

plink --file data/apoly_populations.plink --make-bed --noweb --out data/apoly_popualtions
