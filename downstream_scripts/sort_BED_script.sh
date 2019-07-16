#!/bin/bash
#SBATCH -J sortBED                             # job name
#SBATCH -p super                             # select partion from 128GB, 256GB, 384GB, GPU and super
#SBATCH -N 1                                     # number of nodes requested by user
#SBATCH -t 0-24:00:00                             # run time, format: D-H:M:S (max wallclock time)
#SBATCH -o sortBED.%j.out                         # standard output file name
#SBATCH -e sortBED.%j.err                         # standard error output file name
#SBATCH --mail-user=aishwarya.gogate@utsouthwestern.edu  # specify an email address
#SBATCH --mail-type=BEGIN                               # send email when job status change (start, end, abortion and etc.)
#SBATCH --mail-type=END                               # send email when job status change (start, end, abortion and etc.)
#SBATCH --mail-type=FAIL                               # send email when job status change (start, end, abortion and etc.)


##Usage: sort -k1,1 -k2,2n file.bed > file.sorted.bed

sort -k1,1 -k2,2n ETN.bed > ETN_sorted.bed
