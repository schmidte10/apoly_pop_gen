#!/bin/bash
#PBS -j oe
#PBS -N filtering_snps
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

cd ./2024analysis/gatk

module load vcftools
module load bcftools 

bcftools query -f '%POS\n' gatk.vcf.gz | wc -l 
bcftools query -l gatk.vcf.gz| wc -l 

bcftools filter -g 5 --threads 20 -O z -o gatk.indel5bp.vcf.gz gatk.vcf.gz 

bcftools query -f '%POS\n' gatk.indel5bp.vcf.gz | wc -l 
bcftools query -l gatk.indel5bp.vcf.gz | wc -l 

bcftools filter -e 'QD<2 | QUAL<30 | SOR>3 | FS>60 | MQ<40 | MQRankSum<-12.5 | ReadPosRankSum<-8' gatk.indel5bp.vcf.gz > gatk.indel5b.hardfilter.vcf.gz 

bcftools query -f '%POS\n' gatk.indel5b.hardfilter.vcf.gz | wc -l 
bcftools query -l gatk.indel5b.hardfilter.vcf.gz | wc -l  

vcftools --gzvcf gatk.indel5bp.hardfilter.vcf.gz --relatedness2 --out gatk.indel5bp.hardfilter   

bcftools filter -e 'AC=0 || AC==AN' --threads 20 gatk.indel5bp.hardfilter.vcf.gz |gzip > gatk.filtered.vcf.gz

bcftools query -f '%POS\n' gatk.filtered.vcf.gz | wc -l 
bcftools query -l gatk.filtered.vcf.gz | wc -l  

vcftools --gzvcf gatk.filtered.vcf.gz  --depth 
vcftools --gzvcf gatk.filtered.vcf.gz  --missing-indv 
vcftools --gzvcf gatk.filtered.vcf.gz  --SNPdensity 10000  

#bcftools query -f '%QD\n' gatk.vcf.gz > QD_values.txt
#bcftools query -f '%QUAL\n' gatk.vcf.gz > QUAL_values.txt
#bcftools query -f '%SOR\n' gatk.vcf.gz > SOR_values.txt
#bcftools query -f '%FS\n' gatk.vcf.gz > FS_values.txt
#bcftools query -f '%MQ\n' gatk.vcf.gz > MQ_values.txt
#bcftools query -f '%MQRankSum\n' gatk.vcf.gz > MQRankSum_values.txt
#bcftools query -f '%ReadPosRankSum\n' gatk.vcf.gz > RPRS_values.txt  

bcftools query -f '%QD\n' gatk.filtered.vcf.gz > QD_values_filter.txt
bcftools query -f '%QUAL\n' gatk.filtered.vcf.gz > QUAL_values_filter.txt
bcftools query -f '%SOR\n' gatk.filtered.vcf.gz > SOR_values_filter.txt
bcftools query -f '%FS\n' gatk.filtered.vcf.gz > FS_values_filter.txt
bcftools query -f '%MQ\n' gatk.filtered.vcf.gz > MQ_values_filter.txt
bcftools query -f '%MQRankSum\n' gatk.filtered.vcf.gz > MQRankSum_values_filter.txt
bcftools query -f '%ReadPosRankSum\n' gatk.filtered.vcf.gz > RPRS_values_filter.txt


 

