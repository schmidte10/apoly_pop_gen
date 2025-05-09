#!/bin/bash
#PBS -j oe
#PBS -N PGD.Spider.conv
#PBS -l select=1:ncpus=2:mem=100gb
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

cd ./2024analysis/04.OutlierAnalysis

#java -Xmx100G -Xms10G -jar programs/PGDSpider_2.1.1.5/PGDSpider2-cli.jar -inputfile data/freebayes_f3_chr1.vcf -inputformat VCF -outputfile analysis/bayescan/apoly_populations2.pgd.xml -outputformat PGD -spid analysis/bayescan/nano_PGD.spid
java -Xmx100G -Xms10G -jar programs/PGDSpider_2.1.1.5/PGDSpider2-cli.jar -inputfile analysis/bayescan/apoly_populations2.pgd.xml -inputformat PGD -outputfile analysis/bayescan/apoly_populations2.bs -outputformat GESTE_BAYE_SCAN