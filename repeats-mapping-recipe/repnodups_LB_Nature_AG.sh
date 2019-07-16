#!/bin/bash

#SBATCH --job-name=nodups_mm10repeats
#SBATCH --partition=super
#SBATCH --nodes=1
#SBATCH --time=0-72:00:00
#SBATCH --output=nodups.%j.out
#SBATCH --error=nodups.%j.time
#SBATCH --mail-user=aishwarya.gogate@utsouthwestern.edu
#SBATCH --mail-type=ALL

##LB Nature paper H3K9me3 ChIP data
##Note: This published data is single end. Change -s to ‘paired-end' for PE data.

/home2/s171755/chip-seq-recipe/remove-duplicates.sh -s 'single-end' -a '/project/GCRB/Banaszynski_lab/shared/Yu-Ching/Repeat_Analysis/ChIP_data/LB_Nature/bam_files_AG/align-bowtie-se.sh-1.0.0/ESCH33WTH3K9me3.fastq.gz.sorted.bam' -o '/project/GCRB/Banaszynski_lab/shared/Yu-Ching/Repeat_Analysis/ChIP_data/LB_Nature/nodup_files_AG'
/home2/s171755/chip-seq-recipe/remove-duplicates.sh -s 'single-end' -a '/project/GCRB/Banaszynski_lab/shared/Yu-Ching/Repeat_Analysis/ChIP_data/LB_Nature/bam_files_AG/align-bowtie-se.sh-1.0.0/ESCH33KO1H3K9me3.fastq.gz.sorted.bam' -o '/project/GCRB/Banaszynski_lab/shared/Yu-Ching/Repeat_Analysis/ChIP_data/LB_Nature/nodup_files_AG'
/home2/s171755/chip-seq-recipe/remove-duplicates.sh -s 'single-end' -a '/project/GCRB/Banaszynski_lab/shared/Yu-Ching/Repeat_Analysis/ChIP_data/LB_Nature/bam_files_AG/align-bowtie-se.sh-1.0.0/ESCH33KO2H3K9me3.fastq.gz.sorted.bam' -o '/project/GCRB/Banaszynski_lab/shared/Yu-Ching/Repeat_Analysis/ChIP_data/LB_Nature/nodup_files_AG'
/home2/s171755/chip-seq-recipe/remove-duplicates.sh -s 'single-end' -a '/project/GCRB/Banaszynski_lab/shared/Yu-Ching/Repeat_Analysis/ChIP_data/LB_Nature/bam_files_AG/align-bowtie-se.sh-1.0.0/ESCH33WTInput.fastq.gz.sorted.bam' -o '/project/GCRB/Banaszynski_lab/shared/Yu-Ching/Repeat_Analysis/ChIP_data/LB_Nature/nodup_files_AG'
/home2/s171755/chip-seq-recipe/remove-duplicates.sh -s 'single-end' -a '/project/GCRB/Banaszynski_lab/shared/Yu-Ching/Repeat_Analysis/ChIP_data/LB_Nature/bam_files_AG/align-bowtie-se.sh-1.0.0/ESCH33KO1Input.fastq.gz.sorted.bam' -o '/project/GCRB/Banaszynski_lab/shared/Yu-Ching/Repeat_Analysis/ChIP_data/LB_Nature/nodup_files_AG'
/home2/s171755/chip-seq-recipe/remove-duplicates.sh -s 'single-end' -a '/project/GCRB/Banaszynski_lab/shared/Yu-Ching/Repeat_Analysis/ChIP_data/LB_Nature/bam_files_AG/align-bowtie-se.sh-1.0.0/ESCH33KO2Input.fastq.gz.sorted.bam' -o '/project/GCRB/Banaszynski_lab/shared/Yu-Ching/Repeat_Analysis/ChIP_data/LB_Nature/nodup_files_AG'
/home2/s171755/chip-seq-recipe/remove-duplicates.sh -s 'single-end' -a '/project/GCRB/Banaszynski_lab/shared/Yu-Ching/Repeat_Analysis/ChIP_data/LB_Nature/bam_files_AG/align-bowtie-se.sh-1.0.0/ESCDAXXWTH3K9me3.fastq.gz.sorted.bam' -o '/project/GCRB/Banaszynski_lab/shared/Yu-Ching/Repeat_Analysis/ChIP_data/LB_Nature/nodup_files_AG'
/home2/s171755/chip-seq-recipe/remove-duplicates.sh -s 'single-end' -a '/project/GCRB/Banaszynski_lab/shared/Yu-Ching/Repeat_Analysis/ChIP_data/LB_Nature/bam_files_AG/align-bowtie-se.sh-1.0.0/ESCDAXXnullH3K9me3.fastq.gz.sorted.bam' -o '/project/GCRB/Banaszynski_lab/shared/Yu-Ching/Repeat_Analysis/ChIP_data/LB_Nature/nodup_files_AG'
/home2/s171755/chip-seq-recipe/remove-duplicates.sh -s 'single-end' -a '/project/GCRB/Banaszynski_lab/shared/Yu-Ching/Repeat_Analysis/ChIP_data/LB_Nature/bam_files_AG/align-bowtie-se.sh-1.0.0/ESCATRXWTH3K9me3.fastq.gz.sorted.bam' -o '/project/GCRB/Banaszynski_lab/shared/Yu-Ching/Repeat_Analysis/ChIP_data/LB_Nature/nodup_files_AG'
/home2/s171755/chip-seq-recipe/remove-duplicates.sh -s 'single-end' -a '/project/GCRB/Banaszynski_lab/shared/Yu-Ching/Repeat_Analysis/ChIP_data/LB_Nature/bam_files_AG/align-bowtie-se.sh-1.0.0/ESCATRXnullH3K9me3.fastq.gz.sorted.bam' -o '/project/GCRB/Banaszynski_lab/shared/Yu-Ching/Repeat_Analysis/ChIP_data/LB_Nature/nodup_files_AG'
/home2/s171755/chip-seq-recipe/remove-duplicates.sh -s 'single-end' -a '/project/GCRB/Banaszynski_lab/shared/Yu-Ching/Repeat_Analysis/ChIP_data/LB_Nature/bam_files_AG/align-bowtie-se.sh-1.0.0/ESCDAXXWTInput.fastq.gz.sorted.bam' -o '/project/GCRB/Banaszynski_lab/shared/Yu-Ching/Repeat_Analysis/ChIP_data/LB_Nature/nodup_files_AG'
/home2/s171755/chip-seq-recipe/remove-duplicates.sh -s 'single-end' -a '/project/GCRB/Banaszynski_lab/shared/Yu-Ching/Repeat_Analysis/ChIP_data/LB_Nature/bam_files_AG/align-bowtie-se.sh-1.0.0/ESCDAXXnullInput.fastq.gz.sorted.bam' -o '/project/GCRB/Banaszynski_lab/shared/Yu-Ching/Repeat_Analysis/ChIP_data/LB_Nature/nodup_files_AG'
/home2/s171755/chip-seq-recipe/remove-duplicates.sh -s 'single-end' -a '/project/GCRB/Banaszynski_lab/shared/Yu-Ching/Repeat_Analysis/ChIP_data/LB_Nature/bam_files_AG/align-bowtie-se.sh-1.0.0/ESCATRXWTInput.fastq.gz.sorted.bam' -o '/project/GCRB/Banaszynski_lab/shared/Yu-Ching/Repeat_Analysis/ChIP_data/LB_Nature/nodup_files_AG'
/home2/s171755/chip-seq-recipe/remove-duplicates.sh -s 'single-end' -a '/project/GCRB/Banaszynski_lab/shared/Yu-Ching/Repeat_Analysis/ChIP_data/LB_Nature/bam_files_AG/align-bowtie-se.sh-1.0.0/ESCATRXnullInput.fastq.gz.sorted.bam' -o '/project/GCRB/Banaszynski_lab/shared/Yu-Ching/Repeat_Analysis/ChIP_data/LB_Nature/nodup_files_AG'
