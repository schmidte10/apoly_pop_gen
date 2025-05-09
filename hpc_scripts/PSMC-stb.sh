#!/bin/bash
#PBS -j oe 
#PBS -J 1-4
#PBS -N PSMC-stb
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

bam_files=(10.PSMC/stbees/*_marked_sorted.bam)
bam=${bam_files[$PBS_ARRAY_INDEX - 1]}
sample=$(basename "$bam" _marked_sorted.bam) 

for n in {1..10}; do
	psmc -b -p 4+25*2+4+6 10.PSMC/stbees/${sample}.split.psmcfa -o 10.PSMC/stbees/bootstrap/${sample}.round-${n}.psmc  
done

if [[ $PBS_ARRAY_INDEX -eq 8 ]]; then
cat 10.PSMC/stbees/bootstrap/${sample}.round-*.psmc > 10.PSMC/stbees/${sample}_bs_combined.psmc 
rm 10.PSMC/stbees/bootstrap/${sample}.round-*.psmc
fi
