#!/bin/bash
#PBS -j oe
#PBS -N missing_site
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

module load bedtools2 
module load bioawk

cd 2024analysis 

path=mapped_bams 
declare -A BAM 

BAM=mapped_bams/CSUD006_S1_marked.bam  

#bioawk -c fastx '{print $name"\t1\t"length($seq)}' ../Acanthochromis_polyacanthus_files/reference2024.fna > reference.bed 
#bedtools makewindows -b ../Acanthochromis_polyacanthus_files/reference.bed -w 9999 -s 2000 > ../Acanthochromis_polyacanthus_files/windows_10k.bed

#cal_missing(){
#bedtools genomecov -ibam $BAM -d | awk '{if($3<3) print$0}' | awk '{print $1"\t"$2"\t"$2}' | bedtools merge | bedtools intersect -a ../Acanthochromis_polyacanthus_files/windows_10k.bed -b - -wao | cut -f1,2,3,7 | awk '{sum[$1"\t"$2"\t"$3]+=$4}END{for(i in sum)print i"\t"sum[i]/10000}' | sort -k1,1 -k2,2n -k3,3n > 05.2.Population_structure/sud.windows.missing.txt
#} 

#export -f cal_missing 

#cal_missing sudbury 

bedtools genomecov -ibam mapped_bams/CTON061_S10_marked.bam -d | awk '{if($3<3) print$0}' | awk '{print $1"\t"$2"\t"$2}' | bedtools merge | bedtools intersect -a ../Acanthochromis_polyacanthus_files/windows_10k.bed -b - -wao | cut -f1,2,3,7 | awk '{sum[$1"\t"$2"\t"$3]+=$4}END{for(i in sum)print i"\t"sum[i]/10000}' | sort -k1,1 -k2,2n -k3,3n > 05.2.Population_structure/ton.windows.missing.txt
bedtools genomecov -ibam mapped_bams/CVLA044_S23_marked.bam -d | awk '{if($3<3) print$0}' | awk '{print $1"\t"$2"\t"$2}' | bedtools merge | bedtools intersect -a ../Acanthochromis_polyacanthus_files/windows_10k.bed -b - -wao | cut -f1,2,3,7 | awk '{sum[$1"\t"$2"\t"$3]+=$4}END{for(i in sum)print i"\t"sum[i]/10000}' | sort -k1,1 -k2,2n -k3,3n > 05.2.Population_structure/vla.windows.missing.txt
bedtools genomecov -ibam mapped_bams/LCHA114_S22_marked.bam -d | awk '{if($3<3) print$0}' | awk '{print $1"\t"$2"\t"$2}' | bedtools merge | bedtools intersect -a ../Acanthochromis_polyacanthus_files/windows_10k.bed -b - -wao | cut -f1,2,3,7 | awk '{sum[$1"\t"$2"\t"$3]+=$4}END{for(i in sum)print i"\t"sum[i]/10000}' | sort -k1,1 -k2,2n -k3,3n > 05.2.Population_structure/cha.windows.missing.txt
bedtools genomecov -ibam mapped_bams/LCKM154_S27_marked.bam -d | awk '{if($3<3) print$0}' | awk '{print $1"\t"$2"\t"$2}' | bedtools merge | bedtools intersect -a ../Acanthochromis_polyacanthus_files/windows_10k.bed -b - -wao | cut -f1,2,3,7 | awk '{sum[$1"\t"$2"\t"$3]+=$4}END{for(i in sum)print i"\t"sum[i]/10000}' | sort -k1,1 -k2,2n -k3,3n > 05.2.Population_structure/ckm.windows.missing.txt
bedtools genomecov -ibam mapped_bams/LHER700_S22_marked.bam -d | awk '{if($3<3) print$0}' | awk '{print $1"\t"$2"\t"$2}' | bedtools merge | bedtools intersect -a ../Acanthochromis_polyacanthus_files/windows_10k.bed -b - -wao | cut -f1,2,3,7 | awk '{sum[$1"\t"$2"\t"$3]+=$4}END{for(i in sum)print i"\t"sum[i]/10000}' | sort -k1,1 -k2,2n -k3,3n > 05.2.Population_structure/her.windows.missing.txt
bedtools genomecov -ibam mapped_bams/LKES142_S35_marked.bam -d | awk '{if($3<3) print$0}' | awk '{print $1"\t"$2"\t"$2}' | bedtools merge | bedtools intersect -a ../Acanthochromis_polyacanthus_files/windows_10k.bed -b - -wao | cut -f1,2,3,7 | awk '{sum[$1"\t"$2"\t"$3]+=$4}END{for(i in sum)print i"\t"sum[i]/10000}' | sort -k1,1 -k2,2n -k3,3n > 05.2.Population_structure/kes.windows.missing.txt
bedtools genomecov -ibam mapped_bams/LSBI021_S1_marked.bam -d | awk '{if($3<3) print$0}' | awk '{print $1"\t"$2"\t"$2}' | bedtools merge | bedtools intersect -a ../Acanthochromis_polyacanthus_files/windows_10k.bed -b - -wao | cut -f1,2,3,7 | awk '{sum[$1"\t"$2"\t"$3]+=$4}END{for(i in sum)print i"\t"sum[i]/10000}' | sort -k1,1 -k2,2n -k3,3n > 05.2.Population_structure/sbi.windows.missing.txt
bedtools genomecov -ibam mapped_bams/TDUG708_S30_marked.bam -d | awk '{if($3<3) print$0}' | awk '{print $1"\t"$2"\t"$2}' | bedtools merge | bedtools intersect -a ../Acanthochromis_polyacanthus_files/windows_10k.bed -b - -wao | cut -f1,2,3,7 | awk '{sum[$1"\t"$2"\t"$3]+=$4}END{for(i in sum)print i"\t"sum[i]/10000}' | sort -k1,1 -k2,2n -k3,3n > 05.2.Population_structure/tor.windows.missing.txt
bedtools genomecov -ibam mapped_bams/XRHQ997_S13_marked.bam -d | awk '{if($3<3) print$0}' | awk '{print $1"\t"$2"\t"$2}' | bedtools merge | bedtools intersect -a ../Acanthochromis_polyacanthus_files/windows_10k.bed -b - -wao | cut -f1,2,3,7 | awk '{sum[$1"\t"$2"\t"$3]+=$4}END{for(i in sum)print i"\t"sum[i]/10000}' | sort -k1,1 -k2,2n -k3,3n > 05.2.Population_structure/rhq.windows.missing.txt
