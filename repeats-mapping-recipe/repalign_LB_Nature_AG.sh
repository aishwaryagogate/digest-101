#!/bin/bash
#SBATCH --job-name=align_mm10repeats                              # job name
#SBATCH --partition=super                                 # select partion from 128GB, 256GB, 384GB, GPU and super
#SBATCH --nodes=1                                         # number of nodes requested by user
#SBATCH --time=0-96:00:00                                 # run time, format: D-H:M:S (max wallclock time)
#SBATCH --output=align.%j.out                         # standard output file name
#SBATCH --error=align.%j.err                         # standard error output file name
#SBATCH --mail-user=aishwarya.gogate@utsouthwestern.edu           # specify an email address
#SBATCH --mail-type=ALL                                   # send email when job status change (start, end, abortion and etc.)

##LB Nature paper H3K9me3 ChIP data
##Made separate bowtie index for single end alignment!!
##Note:This published data is single end. Change script to align-bowtie-pe-mod.sh for paired-end data.


/home2/s171755/chip-seq-recipe/align-bowtie-se-mod.sh -f '/project/GCRB/Banaszynski_lab/shared/Yu-Ching/Repeat_Analysis/ChIP_data/LB_Nature/fastq_files/ESCH33WTH3K9me3.fastq.gz' -r 'mm10repeatsAG' -o '/project/GCRB/Banaszynski_lab/shared/Yu-Ching/Repeat_Analysis/ChIP_data/LB_Nature/bam_files_AG'
/home2/s171755/chip-seq-recipe/align-bowtie-se-mod.sh -f '/project/GCRB/Banaszynski_lab/shared/Yu-Ching/Repeat_Analysis/ChIP_data/LB_Nature/fastq_files/ESCH33KO1H3K9me3.fastq.gz' -r 'mm10repeatsAG' -o '/project/GCRB/Banaszynski_lab/shared/Yu-Ching/Repeat_Analysis/ChIP_data/LB_Nature/bam_files_AG'
/home2/s171755/chip-seq-recipe/align-bowtie-se-mod.sh -f '/project/GCRB/Banaszynski_lab/shared/Yu-Ching/Repeat_Analysis/ChIP_data/LB_Nature/fastq_files/ESCH33KO2H3K9me3.fastq.gz' -r 'mm10repeatsAG' -o '/project/GCRB/Banaszynski_lab/shared/Yu-Ching/Repeat_Analysis/ChIP_data/LB_Nature/bam_files_AG'
/home2/s171755/chip-seq-recipe/align-bowtie-se-mod.sh -f '/project/GCRB/Banaszynski_lab/shared/Yu-Ching/Repeat_Analysis/ChIP_data/LB_Nature/fastq_files/ESCH33WTInput.fastq.gz' -r 'mm10repeatsAG' -o '/project/GCRB/Banaszynski_lab/shared/Yu-Ching/Repeat_Analysis/ChIP_data/LB_Nature/bam_files_AG'
/home2/s171755/chip-seq-recipe/align-bowtie-se-mod.sh -f '/project/GCRB/Banaszynski_lab/shared/Yu-Ching/Repeat_Analysis/ChIP_data/LB_Nature/fastq_files/ESCH33KO1Input.fastq.gz' -r 'mm10repeatsAG' -o '/project/GCRB/Banaszynski_lab/shared/Yu-Ching/Repeat_Analysis/ChIP_data/LB_Nature/bam_files_AG'
/home2/s171755/chip-seq-recipe/align-bowtie-se-mod.sh -f '/project/GCRB/Banaszynski_lab/shared/Yu-Ching/Repeat_Analysis/ChIP_data/LB_Nature/fastq_files/ESCH33KO2Input.fastq.gz' -r 'mm10repeatsAG' -o '/project/GCRB/Banaszynski_lab/shared/Yu-Ching/Repeat_Analysis/ChIP_data/LB_Nature/bam_files_AG'
/home2/s171755/chip-seq-recipe/align-bowtie-se-mod.sh -f '/project/GCRB/Banaszynski_lab/shared/Yu-Ching/Repeat_Analysis/ChIP_data/LB_Nature/fastq_files/ESCDAXXWTH3K9me3.fastq.gz' -r 'mm10repeatsAG' -o '/project/GCRB/Banaszynski_lab/shared/Yu-Ching/Repeat_Analysis/ChIP_data/LB_Nature/bam_files_AG'
/home2/s171755/chip-seq-recipe/align-bowtie-se-mod.sh -f '/project/GCRB/Banaszynski_lab/shared/Yu-Ching/Repeat_Analysis/ChIP_data/LB_Nature/fastq_files/ESCDAXXnullH3K9me3.fastq.gz' -r 'mm10repeatsAG' -o '/project/GCRB/Banaszynski_lab/shared/Yu-Ching/Repeat_Analysis/ChIP_data/LB_Nature/bam_files_AG'
/home2/s171755/chip-seq-recipe/align-bowtie-se-mod.sh -f '/project/GCRB/Banaszynski_lab/shared/Yu-Ching/Repeat_Analysis/ChIP_data/LB_Nature/fastq_files/ESCATRXWTH3K9me3.fastq.gz' -r 'mm10repeatsAG' -o '/project/GCRB/Banaszynski_lab/shared/Yu-Ching/Repeat_Analysis/ChIP_data/LB_Nature/bam_files_AG'
/home2/s171755/chip-seq-recipe/align-bowtie-se-mod.sh -f '/project/GCRB/Banaszynski_lab/shared/Yu-Ching/Repeat_Analysis/ChIP_data/LB_Nature/fastq_files/ESCATRXnullH3K9me3.fastq.gz' -r 'mm10repeatsAG' -o '/project/GCRB/Banaszynski_lab/shared/Yu-Ching/Repeat_Analysis/ChIP_data/LB_Nature/bam_files_AG'
/home2/s171755/chip-seq-recipe/align-bowtie-se-mod.sh -f '/project/GCRB/Banaszynski_lab/shared/Yu-Ching/Repeat_Analysis/ChIP_data/LB_Nature/fastq_files/ESCDAXXWTInput.fastq.gz' -r 'mm10repeatsAG' -o '/project/GCRB/Banaszynski_lab/shared/Yu-Ching/Repeat_Analysis/ChIP_data/LB_Nature/bam_files_AG'
/home2/s171755/chip-seq-recipe/align-bowtie-se-mod.sh -f '/project/GCRB/Banaszynski_lab/shared/Yu-Ching/Repeat_Analysis/ChIP_data/LB_Nature/fastq_files/ESCDAXXnullInput.fastq.gz' -r 'mm10repeatsAG' -o '/project/GCRB/Banaszynski_lab/shared/Yu-Ching/Repeat_Analysis/ChIP_data/LB_Nature/bam_files_AG'
/home2/s171755/chip-seq-recipe/align-bowtie-se-mod.sh -f '/project/GCRB/Banaszynski_lab/shared/Yu-Ching/Repeat_Analysis/ChIP_data/LB_Nature/fastq_files/ESCATRXWTInput.fastq.gz' -r 'mm10repeatsAG' -o '/project/GCRB/Banaszynski_lab/shared/Yu-Ching/Repeat_Analysis/ChIP_data/LB_Nature/bam_files_AG'
/home2/s171755/chip-seq-recipe/align-bowtie-se-mod.sh -f '/project/GCRB/Banaszynski_lab/shared/Yu-Ching/Repeat_Analysis/ChIP_data/LB_Nature/fastq_files/ESCATRXnullInput.fastq.gz' -r 'mm10repeatsAG' -o '/project/GCRB/Banaszynski_lab/shared/Yu-Ching/Repeat_Analysis/ChIP_data/LB_Nature/bam_files_AG'
