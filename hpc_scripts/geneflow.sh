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

cd ./2024analysis

java -Xmx1024m -Xms512m -jar /opt/pgdspider/PGDSpider2-cli.jar -inputfile data/thinned_output.vcf -inputformat VCF -outputfile 09.geneflow/apoly_pop_original.pgd -outputformat PGD -spid 09.geneflow/nano_PGD2.spid
java -Xmx1024m -Xms512m -jar /opt/pgdspider/PGDSpider2-cli.jar -inputfile data/thinned_output.vcf -inputformat VCF -outputfile 09.geneflow/apoly_pop_GENEPOP.txt -outputformat GENEPOP -spid 09.geneflow/VCF_GENEPOP.spid
