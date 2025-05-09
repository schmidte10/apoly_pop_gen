#!/bin/bash
#PBS -j oe 
#PBS -J 1-10
#PBS -N PSMC-vla
#PBS -l select=1:ncpus=2:mem=25gb 
#PBS -l walltime=36:00:00
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
module load samtools
module load psmc

#while read -r old_name new_name; do
  #echo "s/^>$old_name/>$new_name/"
#done < Acanthochromis_polyacanthus_files/chr_names_refmap.txt | sed -f - Acanthochromis_polyacanthus_files/reference2024.fna > Acanthochromis_polyacanthus_files/updated_reference2024.fna

#samtools faidx Acanthochromis_polyacanthus_files/updated_reference2024.fna 

## step below takes ~28 hours on 4 ncpus and 200 gb of memory

#for f in mapped_bams/*_marked.bam;do 
#echo $f; 
#samtools sort $f -o ${f%.bam}_sorted.bam  
#done 

#for f in 10.PSMC/*_marked_sorted.bam;do 
#echo $f
#samtools index  ${f%.bam}.bam 
#done
 


# Code below with 200gb mem and 8 ncore took 72 hours to run ~45-50 samples. - likely need 5 days to run through 83 samples

#for bam in 10.PSMC/test/*_marked_sorted.bam; do
    #sample=$(basename "$bam" _marked_sorted.bam)
    #while read chr; do
    #bcftools mpileup -Q 30 -q 30 -C 50 -f Acanthochromis_polyacanthus_files/updated_reference2024.fna -r $chr 10.PSMC/test/${sample}_marked_sorted.bam | \
    #bcftools call -c | \
    #programs/vcfutils.pl vcf2fq -d 10 -D 80 -Q 30 > 10.PSMC/test/${sample}_${chr}.fq
#done < data/chr_names.txt
#done

#for bam in 10.PSMC/*_marked_sorted.bam; do
    #sample=$(basename "$bam" _marked_sorted.bam) 
    	#cat 10.PSMC/${sample}_*.fq > 10.PSMC/${sample}_consensus.fq 
    	#rm 10.PSMC/${sample}_NC*.fq 
    	#fq2psmcfa 10.PSMC/${sample}_consensus.fq > 10.PSMC/${sample}.psmcfa 
    	#psmc -p 4+25*2+4+6 10.PSMC/${sample}.psmcfa -o 10.PSMC/${sample}.psmc 
    	#splitfa 10.PSMC/${sample}.psmcfa > 10.PSMC/${sample}.split.psmcfa
#done

## SUDBURY
#bam_files=(10.PSMC/*_marked_sorted.bam)
#bam=${bam_files[$PBS_ARRAY_INDEX - 1]}
#sample=$(basename "$bam" _marked_sorted.bam) 

#for n in {1..10}; do
#psmc -b -p 4+25*2+4+6 10.PSMC/${sample}.split.psmcfa -o 10.PSMC/bootstrap/${sample}.round-${n}.psmc  
#done

#if [[ $PBS_ARRAY_INDEX -eq 8 ]]; then
#cat 10.PSMC/bootstrap/${sample}.round-*.psmc > 10.PSMC/${sample}_bs_combined.psmc 
#rm 10.PSMC/bootstrap/${sample}.round-*.psmc
#fi

## VLASSOF 

bam_files=(10.PSMC/vlassof/*_marked_sorted.bam)
bam=${bam_files[$PBS_ARRAY_INDEX - 1]}
sample=$(basename "$bam" _marked_sorted.bam) 

for n in {1..10}; do
	psmc -b -p 4+25*2+4+6 10.PSMC/vlassof/${sample}.split.psmcfa -o 10.PSMC/vlassof/bootstrap/${sample}.round-${n}.psmc  
done

if [[ $PBS_ARRAY_INDEX -eq 8 ]]; then
cat 10.PSMC/vlassof/bootstrap/${sample}.round-*.psmc > 10.PSMC/vlassof/${sample}_bs_combined.psmc 
rm 10.PSMC/vlassof/bootstrap/${sample}.round-*.psmc
fi

#plotting 
#psmc_plot.pl -u 2.58e-8 -g 3 10.PSMC/psmc_plot 10.PSMC/CSUD006_S1_bs_combined.psmc
#ps2pdf 10.PSMC/psmc_plot.eps 10.PSMC/psmc_plot.pdf



#calculate divergence time to per generation mutation rate 
#using rough calculation of 1.1 * 10-7, which is conservative because it is based on hyper-variable regions witin the mitochondrial DNA # van Herwerden and Doherty 2006

#seq 1 50 | parallel -j 2 psmc -b -p 4+25*2+4+6 10.PSMC/CSUD008_S2.split.psmcfa -o 10.PSMC/bootstrap/round-{}.psmc 
