#!/bin/bash
#PBS -j oe
#PBS -N baypass
#PBS -l select=1:ncpus=4:mem=100gb
#PBS -l walltime=200:00:00
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

module load plink
module load baypass_public 

cd ./2025wgs

#awk '{print $3, "\t",$1}' data/poppop.txt > data/apoly_populations_metadata_POPIND.txt

# remove first 2 columns 
#cut -f 3- data/apoly_populations.plink.ped > x.delete 

#paste data/apoly_populations_metadata_POPIND.txt x.delete > analysis/baypass/apoly_populations.plink.ped 
#rm x.delete 

# tr -d '\r' < analysis/baypass/apoly_populations.plink.ped > analysis/baypass/apoly_populations.plink.ped

#plink --file analysis/baypass/apoly_populations.plink --allow-extra-chr --freq counts --family --out analysis/baypass/apoly_populations 

#npop2=20 # number of pop times 2
#tail -n +2 analysis/baypass/apoly_populations.frq.strat | awk '{ $9 = $8 - $7 } 1' | awk '{print $7,$9}' | tr " " "\n" | awk -v pp=${npop2} '{if (NR % pp == 0){a=a $0"";print a; a=""} else a=a $0" "}' > analysis/baypass/apoly_populations_baypass.txt 

#g_baypass singularity run /fast/tmp/containers/baypass_public-v2.41.sif g_baypass -npop 6 -gfile baypass/apoly_populations_baypass.txt -outprefix baypass/apoly_populations_baypass -nthreads 4

g_baypass singularity run /fast/tmp/containers/baypass_public-v2.41.sif g_baypass -npop 6 -gfile baypass/apoly_populations_baypass.txt -efile baypass/metadata/baypass_env_metadata.txt -scalecov -omegafile baypass/apoly_populations_baypass_mat_omega.out -outprefix baypass/apoly_populations_baypass_env -nthreads 4
