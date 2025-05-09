#!/bin/bash
#PBS -j oe
#PBS -N pixy_fst
#PBS -l select=1:ncpus=4:mem=100gb
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

module load pixy

cd ./2024analysis/

pixy --stats fst pi dxy --vcf data/gatk.filtered.relaxed_studywide.vcf.gz --populations data/apoly_populations_metadata_INDREG.txt --window_size 10000 --n_cores 8 --bypass_invariant_check yes
