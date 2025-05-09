#!/bin/bash
#PBS -j oe
#PBS -N PGD.Spider.conv
#PBS -l select=1:ncpus=8:mem=250gb
#PBS -l walltime=50:00:00
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

cd ./2024analysis/04.OutlierAnalysis

gunzip data/gatk.filtered.relaxed_studywide.vcf.gz > data/gatk.filtered.relaxed_studywide.vcf
java -Xmx100G -Xms10G -jar programs/PGDSpider_2.1.1.5/PGDSpider2-cli.jar -inputfile data/gatk.filtered.relaxed_studywide.vcf -inputformat VCF -outputfile analysis/bayescan/apoly_pop_original.pgd.xml -outputformat PGD -spid analysis/bayescan/nano_PGD.spid
java -Xmx100G -Xms10G -jar programs/PGDSpider_2.1.1.5/PGDSpider2-cli.jar -inputfile analysis/bayescan/apoly_pop_original.pgd.xml -inputformat PGD -outputfile analysis/bayescan/apoly_pop_original.bs -outputformat GESTE_BAYE_SCAN