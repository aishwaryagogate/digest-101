#!/bin/bash

#SBATCH --job-name=bed_intersect
#SBATCH --partition=super
#SBATCH --nodes=1
#SBATCH --time=0-24:00:00
#SBATCH --output=bedintersect.%j.out
#SBATCH --error=bedintersect.%j.time
#SBATCH --mail-user=aishwarya.gogate@utsouthwestern.edu
#SBATCH --mail-type=ALL

module load bedtools

#To get common peaks
bedtools intersect -a fileA.bed -b fileB.bed > common_peaks.bed

#To get peaks in fileA only
bedtools intersect -v -a fileA.bed -b fileB.bed > fileApeaksonly.bed

#To get peaks in fileB only
bedtools intersect -v -a fileB.bed -b fileA.bed > fileBpeaksonly.bed
