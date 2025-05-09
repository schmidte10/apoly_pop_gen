#!/bin/bash
#PBS -j oe
#PBS -N mitogenome_map
#PBS -l select=1:ncpus=4:mem=100gb 
#PBS -l walltime=10:00:00
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

cd 2024analysis/mapped_bams/


module load samtools
module load bwa 
module load bcftools 
module load bioawk
module load tabix

#ncpus should be set to 18 when running all

#for input_bam in *_marked.bam;do 
#sample=${input_bam%_marked.bam} 
#samtools fastq -F 1024 ${input_bam} > ${sample}_marked_discarded_singeltons.bam | wc -l 
#samtools collate -O ${sample}_marked_discarded_singeltons.bam | samtools fastq -1 ${sample}_pair_1.fq -2 ${sample}_pair2.fq -0 /dev/null -s /dev/null -n #correct pair2.fq name should include '_'
#bwa mem -t 16 ../mitogenome/Organelle_final_assembly.fa ${sample}_pair_1.fq ${sample}_pair2.fq | samtools view -b -F 4 - > ${sample}_mito.bam #correct for updated pair2.fq name 
#rm ${sample}_marked_discarded_singeltons.bam 
#rm ${sample}_pair_1.fq 
#rm ${sample}_pair2.fq
#done

#for f in *_mito.bam;do 
#echo $f; 
#samtools sort $f -o ${f%.bam}_sorted.bam
#done 

#for f in *_mito_sorted.bam;do 
#mean_cov=$(samtools depth ${f} | awk '{sum+=$3}END{print sum/NR}'); 
#echo $f ${mean_cov}; 
#done > all_covsummary.txt 

#cat all_covsummary.txt | awk '$2>200{print $1}' | sed 's/_mito_sorted.bam//' > goodcov_samples_thres200.txt
#cat all_covsummary.txt | awk '$2>500{print $1}' | sed 's/_mito_sorted.bam//' > goodcov_samples_thres500.txt

call_cons_jia(){
	mitogenome=../mitogenome/Organelle_final_assembly.fa                   #Change to the name of your mitogenome

	bamfile=$1

	sample=${bamfile%_mito_sorted.bam}

	#mean_cov=$(samtools depth ${bamfile}_mito_sorted.bam | awk '{sum += $3}END{print sum/NR}')
	#mean_cov=$(printf "%.0f" $mean_cov)
	#let max_dp=${mean_cov}*2

	#bcftools mpileup -Ou -f ${mitogenome} ${bamfile}_mito_sorted.bam | 
	#bcftools call -c --ploidy 1 - | 
	#bcftools norm -Ob -f ${mitogenome} - | 
	#bcftools filter -Oz --SnpGap 3 --IndelGap 5 -e "QUAL<10 || DP > ${max_dp} || DP<3" - > ${sample}.vcf.gz 

gunzip ${sample}.vcf.gz  
bgzip ${sample}.vcf 
tabix ${sample}.vcf.gz

cat ${mitogenome} | bcftools consensus -a - -e 'type="indel"' ${sample}.vcf.gz |\
	bioawk -c fastx -v samp=${sample} '{printf(">%s\n%s\n", samp, $seq)}' > ${sample}_consensus.fa


}

for goodbam in $(cat goodcov_samples_thres200.txt);do 
call_cons_jia $goodbam

done

#test 
