#!/bin/bash

#SBATCH --job-name rpkm_boxplot
#SBATCH --partition=super
#SBATCH --nodes=2
#SBATCH -t 0-24:0:0
#SBATCH -o rpkm_job_%j.out
#SBATCH -e rpkm_job_%j.err
#SBATCH --mail-type ALL
#SBATCH --mail-user aishwarya.gogate@UTSouthwestern.edu

module load R/3.2.1-intel

ulimit -l unlimited
R --vanilla < boxplot_H33.r

# END OF SCRIPT
