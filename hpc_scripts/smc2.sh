#!/bin/bash
#PBS -j oe 
#PBS -N SMC++
#PBS -l select=1:ncpus=10:mem=100gb 
#PBS -l walltime=70:00:00
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

module load smcpp

#smc++ estimate --cores 10 -o 11.SMC++/smc_files_csud --base SUD --em-iterations 20 --timepoints 20 20000 --thinning 2000 --knots 40 2.85e-8 11.SMC++/smc_files/merged/CSUD.smc.gz 2> error_log.txt 
#smc++ estimate --cores 10 -o 11.SMC++/smc_files_lher --base HER --em-iterations 20 --timepoints 20 20000 --thinning 2000 --knots 40 2.85e-8 11.SMC++/smc_files/merged/LHER.smc.gz 2> error_log_her.txt 
smc++ estimate --cores 10 -o 11.SMC++/smc_files_lckc --base CKC --em-iterations 20 --timepoints 20 20000 --thinning 2000 --knots 40 2.85e-8 11.SMC++/smc_files/merged/LCKC.smc.gz 2> error_log_ckc.txt 


#smc++ plot lher_plot.pdf 11.SMC++/smc_files_lckc/LCKC.final.json
