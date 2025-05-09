#!/bin/bash
#PBS -j oe
#PBS -N fastsimcoal
#PBS -l select=1:ncpus=4:mem=100gb 
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

cd 2024analysis/


module load bedtools2 
module load bcftools 
module load tabix 
module load python3

#awk '{print $1"\t"($4-1)"\t"$5}' genomic.bed > genomic_formatted.bed
#bedtools subtract -a data/gatk.filtered.relaxed_studywide.vcf.gz -b data/genomic_formatted.bed | bgzip > data/gatk.filtered.relaxed_studywide_nongenic.vcf.gz  
#bcftools view -h data/gatk.filtered.relaxed_studywide.vcf.gz > header.txt
#(cat header.txt; zcat gatk.filtered.relaxed_studywide_nongenic.vcf.gz) | bgzip > nongenic.vcf.gz
bcftools +prune -m 0.3 -e 'F_MISSING>=0.01' -w 1000 data/nongenic.vcf.gz -o data/gatk.filtered.relaxed_studywide.r03nomiss.vcf.gz

python3 programs/easySFS/easySFS.py -i data/gatk.filtered.relaxed_studywide.r03nomiss.vcf.gz -p data/apoly_populations_metadata_INDREG.txt -a --GQ 20 --dtype int -o sfs --proj=8,8,8,8,1  -f

