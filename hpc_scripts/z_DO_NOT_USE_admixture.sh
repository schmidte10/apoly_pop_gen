#!/bin/bash
#PBS -j oe
#PBS -N admixture2
#PBS -l select=1:ncpus=8:mem=100gb
#PBS -l walltime=100:00:00
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

module load vcftools 
module load plink 
module load admixture

cd ./2024analysis/

#vcftools --gzvcf data/gatk.filtered.relaxed_studywide.vcf.gz --plink-tped --out 05.Admixture/apoly_admixture 
#plink --tped 05.Admixture/apoly_admixture.tped --tfam 05.Admixture/apoly_admixture.tfam --make-bed --out 05.Admixture/apoly_admixture 

#for K in echo $(seq 11) ; do admixture --cv=10 -B200 -j8 05.Admixture/apoly_admixture.bed $K | tee log${K}.out; done

admixture --cv=10 -B100 -j8 05.Admixture/apoly_admixture.ped 2 | tee log2.out
admixture --cv=10 -B100 -j8 05.Admixture/apoly_admixture.ped 3 | tee log3.out
admixture --cv=10 -B100 -j8 05.Admixture/apoly_admixture.ped 4 | tee log4.out
admixture --cv=10 -B100 -j8 05.Admixture/apoly_admixture.ped 5 | tee log5.out
admixture --cv=10 -B100 -j8 05.Admixture/apoly_admixture.ped 6 | tee log6.out
admixture --cv=10 -B100 -j8 05.Admixture/apoly_admixture.ped 7 | tee log7.out
admixture --cv=10 -B100 -j8 05.Admixture/apoly_admixture.ped 8 | tee log8.out
admixture --cv=10 -B100 -j8 05.Admixture/apoly_admixture.ped 9 | tee log9.out
admixture --cv=10 -B100 -j8 05.Admixture/apoly_admixture.ped 10 | tee log10.out

