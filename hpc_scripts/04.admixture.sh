#!/bin/bash
#PBS -j oe
#PBS -N admix_baker
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

module load bcftools 
module load vcftools 
module load plink
module load plink2 
module load admixture

cd ./2024analysis/ 

vcftools --gzvcf gatk/gatk.indel5bp.hardfilter.vcf.gz\
--max-missing 0.5 --min-meanDP 10 --max-meanDP 30\
--minDP 3 --minGQ 40 --mac 2 --remove-filtered-geno-all --recode\
--recode-INFO-all --stdout |\
gzip > 05.2.Population_structure/gatk.filtered.relaxed_studywide.vcf.gz  

vcftools --remove-indv CTON069_S13 --gzvcf 05.2.Population_structure/gatk.filtered.relaxed_studywide.vcf.gz\
--recode --recode-INFO-all --stdout | gzip > 05.2.Population_structure/gatk.filtered.relaxed_studywide_invariants.vcf.gz 

cd 05.2.Population_structure #[change to your own directory]

bcftools view -H gatk.filtered.relaxed_studywide_invariants.vcf.gz | cut -f 1 | uniq | awk '{print $0"\t"$0 }' > relaxed_studywide_invariants.chrom-map.txt

plink2 --vcf ../data/gatk.filtered.relaxed_studywide.vcf.gz --out apolyacanthus --allow-extra-chr --max-alleles 2 --make-bed --set-missing-var-ids @:# --double-id --geno 0.1 --mac 2 --hwe 0.0001 

plink2 --bfile apolyacanthus --indep-pairwise 50 10 0.1 --out apolyacanthus\
    --allow-extra-chr --bad-ld --set-missing-var-ids @:#

plink2 --bfile apolyacanthus --extract apolyacanthus.prune.in --recode ped --allow-extra-chr --out apolyacanthus_ldpruned  

for K in {1..10} 
do
  admixture --cv=10 apolyacanthus_ldpruned.ped $K | tee log.${K}.out
  awk -v K=$K '$1=="CV"{print K, $4}' log.${K}.out >> CV.txt
done








