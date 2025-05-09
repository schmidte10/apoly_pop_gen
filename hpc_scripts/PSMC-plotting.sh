#!/bin/bash
#PBS -j oe 
#PBS -N PSMC-plotting
#PBS -l select=1:ncpus=4:mem=100gb 
#PBS -l walltime=06:00:00
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

module load psmc

## Chauvel
for bam in 10.PSMC/chauvel/*_marked_sorted.bam; do
    sample=$(basename "$bam" _marked_sorted.bam) 
    	#cat 10.PSMC/chauvel/bootstrap/${sample}.round-*.psmc > 10.PSMC/chauvel/${sample}_bs_combined.psmc 
        #rm 10.PSMC/chauvel/bootstrap/${sample}.round-*.psmc 
        psmc_plot.pl -u 2.58e-8 -g 3 -R 10.PSMC/chauvel/${sample} 10.PSMC/chauvel/${sample}_bs_combined.psmc
        ps2pdf 10.PSMC/chauvel/${sample}.eps 10.PSMC/chauvel/${sample}.pdf 
        cat 10.PSMC/chauvel/${sample}.*.txt > 10.PSMC/plot_files/${sample}.txt
        rm 10.PSMC/chauvel/${sample}.*.txt
done

## Heron 
for bam in 10.PSMC/heron/*_marked_sorted.bam; do
    sample=$(basename "$bam" _marked_sorted.bam) 
    	#cat 10.PSMC/heron/bootstrap/${sample}.round-*.psmc > 10.PSMC/heron/${sample}_bs_combined.psmc 
        #rm 10.PSMC/heron/bootstrap/${sample}.round-*.psmc
        psmc_plot.pl -u 2.58e-8 -g 3 -R 10.PSMC/heron/${sample} 10.PSMC/heron/${sample}_bs_combined.psmc
        ps2pdf 10.PSMC/heron/${sample}.eps 10.PSMC/heron/${sample}.pdf
        cat 10.PSMC/heron/${sample}.*.txt > 10.PSMC/plot_files/${sample}.txt
        rm 10.PSMC/heron/${sample}.*.txt
done

## Cockermouth
for bam in 10.PSMC/cockermouth/*_marked_sorted.bam; do
    sample=$(basename "$bam" _marked_sorted.bam) 
    	#cat 10.PSMC/cockermouth/bootstrap/${sample}.round-*.psmc > 10.PSMC/cockermouth/${sample}_bs_combined.psmc 
        #rm 10.PSMC/cockermouth/bootstrap/${sample}.round-*.psmc
        psmc_plot.pl -u 2.58e-8 -g 3 -R 10.PSMC/cockermouth/${sample} 10.PSMC/cockermouth/${sample}_bs_combined.psmc
        ps2pdf 10.PSMC/cockermouth/${sample}.eps 10.PSMC/cockermouth/${sample}.pdf 
        cat 10.PSMC/cockermouth/${sample}.*.txt > 10.PSMC/plot_files/${sample}.txt
        rm 10.PSMC/cockermouth/${sample}.*.txt
done

## St. Bees 
for bam in 10.PSMC/stbees/*_marked_sorted.bam; do
    sample=$(basename "$bam" _marked_sorted.bam) 
    	#cat 10.PSMC/stbees/bootstrap/${sample}.round-*.psmc > 10.PSMC/stbees/${sample}_bs_combined.psmc 
        #rm 10.PSMC/stbees/bootstrap/${sample}.round-*.psmc
        psmc_plot.pl -u 2.58e-8 -g 3 -R 10.PSMC/stbees/${sample} 10.PSMC/stbees/${sample}_bs_combined.psmc
        ps2pdf 10.PSMC/stbees/${sample}.eps 10.PSMC/stbees/${sample}.pdf 
        cat 10.PSMC/stbees/${sample}.*.txt > 10.PSMC/plot_files/${sample}.txt
        rm 10.PSMC/stbees/${sample}.*.txt
done

## Keswick
for bam in 10.PSMC/keswick/*_marked_sorted.bam; do
    sample=$(basename "$bam" _marked_sorted.bam) 
    	#cat 10.PSMC/keswick/bootstrap/${sample}.round-*.psmc > 10.PSMC/keswick/${sample}_bs_combined.psmc 
        #rm 10.PSMC/keswick/bootstrap/${sample}.round-*.psmc
        psmc_plot.pl -u 2.58e-8 -g 3 -R 10.PSMC/keswick/${sample} 10.PSMC/keswick/${sample}_bs_combined.psmc
        ps2pdf 10.PSMC/keswick/${sample}.eps 10.PSMC/keswick/${sample}.pdf 
        cat 10.PSMC/keswick/${sample}.*.txt > 10.PSMC/plot_files/${sample}.txt
        rm 10.PSMC/keswick/${sample}.*.txt
done

## Sudbury
for bam in 10.PSMC/*_marked_sorted.bam; do
    sample=$(basename "$bam" _marked_sorted.bam) 
    	#cat 10.PSMC/bootstrap/${sample}.round-*.psmc > 10.PSMC/${sample}_bs_combined.psmc 
        #rm 10.PSMC/bootstrap/${sample}.round-*.psmc
        psmc_plot.pl -u 2.58e-8 -g 3 -R 10.PSMC/${sample} 10.PSMC/${sample}_bs_combined.psmc
        ps2pdf 10.PSMC/${sample}.eps 10.PSMC/${sample}.pdf
        cat 10.PSMC/${sample}.*.txt > 10.PSMC/plot_files/${sample}.txt
        rm 10.PSMC/${sample}.*.txt
done

## Vlassof
for bam in 10.PSMC/vlassof/*_marked_sorted.bam; do
    sample=$(basename "$bam" _marked_sorted.bam) 
    	#cat 10.PSMC/vlassof/bootstrap/${sample}.round-*.psmc > 10.PSMC/vlassof/${sample}_bs_combined.psmc 
        #rm 10.PSMC/vlassof/bootstrap/${sample}.round-*.psmc
        psmc_plot.pl -u 2.58e-8 -g 3 -R 10.PSMC/vlassof/${sample} 10.PSMC/vlassof/${sample}_bs_combined.psmc
        ps2pdf 10.PSMC/vlassof/${sample}.eps 10.PSMC/vlassof/${sample}.pdf 
        cat 10.PSMC/vlassof/${sample}.*.txt > 10.PSMC/plot_files/${sample}.txt
        rm 10.PSMC/vlassof/${sample}.*.txt
done

## Tongue
for bam in 10.PSMC/tongue/*_marked_sorted.bam; do
    sample=$(basename "$bam" _marked_sorted.bam) 
    	#cat 10.PSMC/tongue/bootstrap/${sample}.round-*.psmc > 10.PSMC/tongue/${sample}_bs_combined.psmc 
        #rm 10.PSMC/tongue/bootstrap/${sample}.round-*.psmc
        psmc_plot.pl -u 2.58e-8 -g 3 -R 10.PSMC/tongue/${sample} 10.PSMC/tongue/${sample}_bs_combined.psmc
        ps2pdf 10.PSMC/tongue/${sample}.eps 10.PSMC/tongue/${sample}.pdf 
        cat 10.PSMC/tongue/${sample}.*.txt > 10.PSMC/plot_files/${sample}.txt
        rm 10.PSMC/tongue/${sample}.*.txt
done

## Torres
for bam in 10.PSMC/torres/*_marked_sorted.bam; do
    sample=$(basename "$bam" _marked_sorted.bam) 
    	#cat 10.PSMC/torres/bootstrap/${sample}.round-*.psmc > 10.PSMC/torres/${sample}_bs_combined.psmc 
        #rm 10.PSMC/torres/bootstrap/${sample}.round-*.psmc
        psmc_plot.pl -u 2.58e-8 -g 3 -R 10.PSMC/torres/${sample} 10.PSMC/torres/${sample}_bs_combined.psmc
        ps2pdf 10.PSMC/torres/${sample}.eps 10.PSMC/torres/${sample}.pdf 
        cat 10.PSMC/torres/${sample}.*.txt > 10.PSMC/plot_files/${sample}.txt
        rm 10.PSMC/torres/${sample}.*.txt
done
