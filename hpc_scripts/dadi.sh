#!/bin/bash
#PBS -j oe
#PBS -N dadi
#PBS -l select=1:ncpus=4:mem=300gb 
#PBS -l walltime=5:00:00
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

cd 2024analysis

singularity run $SING/dadi-2.3.3.sif python3 python_scripts/dadi_script.py

# create sample dictionary 
#ss ={'Sudbury':9, 'Tongue':6, 'Vlassof':8,'Chauvel':8,'Cockermouth':8,'Heron':8,'Keswick':8,'StBees':8,'Torres':15,'ReefHQ':3} 
# pass it as an additional argument 
#dd =dadi.Misc.make_data_dict_vcf("gatk.filetered.relaxed_studywide.vcf.gz", "apoly_populations_metadata_INDPOP.txt", subsample=ss) 
#fs=dadi.Spectrum.from_data_dict(dd, pop_id = ['Sudbury','Tongue','Vlassof','Chauvel','Cockermouth','Heron','Keswick','StBees','Torres','ReefHQ'], projections =[8,6,8,8,8,8,8,8,12,2], polarized = FALSE) 


