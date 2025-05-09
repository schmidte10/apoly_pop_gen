#!/bin/bash
#PBS -j oe
#PBS -N phasing
#PBS -l select=1:ncpus=4:mem=200gb 
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

cd 2024analysis/


module load shapeit5
module load bcftools 

# make .vcf file into .bcf file 
bcftools view -O b -o <output>.bcf <input>.vcf.gz # replace file names

#example code
#phase_common_static singularity run /fast/tmp/containers/shapeit5-5.1.1.sif phase_common_static --input 07.Phasing/gatk.filtered.relaxed_studywide.bcf --filter-maf 0.001 --region NC_067113.1 --output 07.Phasing/phased_NC067113_test.bcf --thread 4

region_file=data/chr_names.txt

while IFS= read -r region 
do 
output_file="07.Phasing/phased_${region}.bcf" 
singularity run /fast/tmp/containers/shapeit5-5.1.1.sif phase_common_static --input data/gatk.filtered.relaxed_studywide.bcf --filter-maf 0.001 --region $region --output $output_file --thread 4 
echo "Finished phasing region $region" 
done < "data/chr_names.txt"


# note for some reason this code missing the last line for chromosome 24 (NC_067136.1).  

ls -1v 07.Phasing/phased_NC_*.bcf > 07.Phasing/tmp/files.txt 
bcftools concat --naive -f 07.Phasing/tmp/files.txt -o 07.Phasing/phased_data/target.phased.bcf --threads 4 
bcftools index -f 07.Phasing/phased_data/target.phased.bcf 

cd 07.Phasing/phased_data

bcftools view -O z -o gatk.filtered.relaxed_studywide.phased.vcf.gz target.phased.bcf