#!/bin/bash
#PBS -j oe 
#PBS -N SMC++
#PBS -l select=1:ncpus=8:mem=100gb 
#PBS -l walltime=12:00:00
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

module load bcftools 
module load tabix
module load bedtools2
module load samtools 
module load smcpp

#while read chr;do 
#bcftools view -r "$chr" -Oz -o 11.SMC++/"${chr}".vcf.gz data/gatk.filtered.relaxed_studywide.vcf.gz 
#tabix 11.SMC++/"${chr}".vcf.gz 
#done < data/chr_names.txt

#for f in mapped_bams/*_marked.bam;do 
#echo $f; 
#samtools sort $f -o ${f%.bam}_sorted.bam  
#done 

#for f in mapped_bams/*_marked_sorted.bam;do 
#echo $f
#samtools index  ${f%.bam}.bam 
#done

#while read chr;do
#samtools depth -r "${chr}" -aa -f data/bamfiles.txt |\
#awk '{sum=0; for(i=3; i<=NF; i++) {sum+=$i}; print $1"\t"$2"\t"sum }'|\
#awk '{if($3<3 || $3>3000) print $1"\t"$2"\t"$2}' |\
#bedtools merge -i stdin | bgzip > 11.SMC++/"${chr}".low_masked.bed.gz
#done < data/chr_names.txt

#zcat *.low_masked.bed.gz | bgzip > all.masked.bed.gz 

#smc++ vcf2smc -d CSUD006_S1 CSUD006_S1\
#--mask 11.SMC++/all.masked.bed.gz data/gatk.filtered.relaxed_studywide.vcf.gz NC_067113.1_CSUD0061_S1.smc.gz NC_067113.1\
#SUD:$(cat data/population_lists/sud.txt | paste -s -d ',') 


#smc++ vcf2smc -d CSUD006_S1 CSUD006_S1\
#--mask 11.SMC++/all.masked.bed.gz data/gatk.filtered.relaxed_studywide.vcf.gz NC_067113.1_CSUD0061_S1.smc.gz NC_067113.1\
#SUD:$(cat data/population_lists/sud.txt | paste -s -d ',') 

chr_file="data/chr_names.txt" 
#sample_file="data/population_lists/sud.txt" 

#for chr in $(cat $chr_file); do
  #for sample in $(cat $sample_file); do
    #smc++ vcf2smc -d $sample $sample \
      #--mask 11.SMC++/all.masked.bed.gz data/gatk.filtered.relaxed_studywide.vcf.gz 11.SMC++/smc_files/${sample}_${chr}.smc.gz $chr \
      #SUD:$(cat $sample_file | paste -s -d ',')
  #done
#done 

#sample_file="data/population_lists/her.txt" 

#for chr in $(cat $chr_file); do
  #for sample in $(cat $sample_file); do
    #smc++ vcf2smc -d $sample $sample \
      #--mask 11.SMC++/all.masked.bed.gz data/gatk.filtered.relaxed_studywide.vcf.gz 11.SMC++/smc_files/her/${sample}_${chr}.smc.gz $chr \
      #HER:$(cat $sample_file | paste -s -d ',')
  #done
#done 

#cat 11.SMC++/smc_files/her/LHER*.smc.gz > 11.SMC++/smc_files/merged/LHER.smc.gz

sample_file="data/population_lists/ckc.txt" 

for chr in $(cat $chr_file); do
  for sample in $(cat $sample_file); do
    smc++ vcf2smc -d $sample $sample \
      --mask 11.SMC++/all.masked.bed.gz data/gatk.filtered.relaxed_studywide.vcf.gz 11.SMC++/smc_files/ckc/${sample}_${chr}.smc.gz $chr \
      CKC:$(cat $sample_file | paste -s -d ',')
  done
done 

cat 11.SMC++/smc_files/ckc/LHER*.smc.gz > 11.SMC++/smc_files/merged/LCKC.smc.gz


