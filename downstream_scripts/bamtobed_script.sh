#!/bin/bash
#SBATCH -J getBED                             # job name
#SBATCH -p super                             # select partion from 128GB, 256GB, 384GB, GPU and super
#SBATCH -N 1                                     # number of nodes requested by user
#SBATCH -t 0-24:00:00                             # run time, format: D-H:M:S (max wallclock time)
#SBATCH -o getBED.%j.out                         # standard output file name
#SBATCH -e getBED.%j.err                         # standard error output file name
#SBATCH --mail-user=aishwarya.gogate@utsouthwestern.edu  # specify an email address
#SBATCH --mail-type=BEGIN                               # send email when job status change (start, end, abortion and etc.)
#SBATCH --mail-type=END                               # send email when job status change (start, end, abortion and etc.)
#SBATCH --mail-type=FAIL                               # send email when job status change (start, end, abortion and etc.)

module load bedtools/2.17.0
module load picard/1.127
module load samtools

##Usage##
bamToBed -i /project/GCRB/Banaszynski_lab/shared/yourfolder/filename.bam > filename.bed

##This outputs an unsorted BED file. You will then have to sort it. Another script is available.
