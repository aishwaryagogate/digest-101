#!/bin/bash
#SBATCH -J RepEnrich_step4                           			 # job name
#SBATCH -p super                                 # select partion from 128GB, 256GB, 384GB, GPU and super
#SBATCH -N 1                                     			 # number of nodes requested by user
#SBATCH -t  0-240:00:00                                 # run time, format: D-H:M:S (max wallclock time)
#SBATCH --output=RepEnrich_step4.%j.out                         # standard output file name
#SBATCH --error=RepEnrich_step4.%j.err                         # standard error output file name
#SBATCH --mail-user=aishwarya.gogate@utsouthwestern.edu 		 # specify an email address
#SBATCH --mail-type=BEGIN                               		 # send email when job status change (start, end, abortion and etc.)
#SBATCH --mail-type=END                               		 # send email when job status change (start, end, abortion and etc.)
#SBATCH --mail-type=FAIL                              		 # send email when job status change (start, end, abortion and etc.)


##Loading modules##
module load iGenomes/2013-03-25
module load bedtools/2.26.0
module load bowtie2/2.2.8-intel
module load samtools/1.6
module load perl/5.28.0

##NOTE: This step takes ages to run, so you might want to run one sample in each script (with 1 node and 240 hrs each).##

##Step 4) Run RepEnrich2
###After we get the unique and multimap fastq, we run RepEnrich2 on the data
###Note : I have modified the original Repenrich2 script for compatibility with BioHPC
python2.7 /work/GCRB/s157114/softwares/RepEnrich2/RepEnrich2_modified.py /work/GCRB/s157114/softwares/repeatmasker_seq/mm10_repeatmasker_without_simple_repeats.txt ./ESC_H33_WT_H3K9me3 ESC_H33_WT_H3K9me3 /work/GCRB/s157114/softwares/repeatmasker_seq/setup_folder_mm10_no_simple_repeats /project/GCRB/Banaszynski_lab/shared/Yu-Ching/Repeat_Analysis/ChIP_data/LB_Nature/RepEnrich_analysis/ESC_H33_WT_H3K9me3_R1_001_multimap.fastq /project/GCRB/Banaszynski_lab/shared/Yu-Ching/Repeat_Analysis/ChIP_data/LB_Nature/RepEnrich_analysis/ESC_H33_WT_H3K9me3_R1_001_unique.bam --pairedend FALSE --allcountmethod TRUE

##Step 4) Run RepEnrich2
###After we get the unique and multimap fastq, we run RepEnrich2 on the data
###Note : I have modified the original Repenrich2 script for compatibility with BioHPC
python2.7 /work/GCRB/s157114/softwares/RepEnrich2/RepEnrich2_modified.py /work/GCRB/s157114/softwares/repeatmasker_seq/mm10_repeatmasker_without_simple_repeats.txt ./ESC_H33_WT_Input ESC_H33_WT_Input /work/GCRB/s157114/softwares/repeatmasker_seq/setup_folder_mm10_no_simple_repeats /project/GCRB/Banaszynski_lab/shared/Yu-Ching/Repeat_Analysis/ChIP_data/LB_Nature/RepEnrich_analysis/ESC_H33_WT_Input_R1_001_multimap.fastq /project/GCRB/Banaszynski_lab/shared/Yu-Ching/Repeat_Analysis/ChIP_data/LB_Nature/RepEnrich_analysis/ESC_H33_WT_Input_R1_001_unique.bam --pairedend FALSE --allcountmethod TRUE

#Step 4) Run RepEnrich2
###After we get the unique and multimap fastq, we run RepEnrich2 on the data
###Note : I have modified the original Repenrich2 script for compatibility with BioHPC
python2.7 /work/GCRB/s157114/softwares/RepEnrich2/RepEnrich2_modified.py /work/GCRB/s157114/softwares/repeatmasker_seq/mm10_repeatmasker_without_simple_repeats.txt ./ESC_H33_KO1_H3K9me3 ESC_H33_KO1_H3K9me3 /work/GCRB/s157114/softwares/repeatmasker_seq/setup_folder_mm10_no_simple_repeats /project/GCRB/Banaszynski_lab/shared/Yu-Ching/Repeat_Analysis/ChIP_data/LB_Nature/RepEnrich_analysis/ESC_H33_KO1_H3K9me3_R1_001_multimap.fastq /project/GCRB/Banaszynski_lab/shared/Yu-Ching/Repeat_Analysis/ChIP_data/LB_Nature/RepEnrich_analysis/ESC_H33_KO1_H3K9me3_R1_001_unique.bam --pairedend FALSE --allcountmethod TRUE

##Step 4) Run RepEnrich2
###After we get the unique and multimap fastq, we run RepEnrich2 on the data
###Note : I have modified the original Repenrich2 script for compatibility with BioHPC
python2.7 /work/GCRB/s157114/softwares/RepEnrich2/RepEnrich2_modified.py /work/GCRB/s157114/softwares/repeatmasker_seq/mm10_repeatmasker_without_simple_repeats.txt ./ESC_H33_KO2_H3K9me3 ESC_H33_KO2_H3K9me3 /work/GCRB/s157114/softwares/repeatmasker_seq/setup_folder_mm10_no_simple_repeats /project/GCRB/Banaszynski_lab/shared/Yu-Ching/Repeat_Analysis/ChIP_data/LB_Nature/RepEnrich_analysis/ESC_H33_KO2_H3K9me3_R1_001_multimap.fastq /project/GCRB/Banaszynski_lab/shared/Yu-Ching/Repeat_Analysis/ChIP_data/LB_Nature/RepEnrich_analysis/ESC_H33_KO2_H3K9me3_R1_001_unique.bam --pairedend FALSE --allcountmethod TRUE

##Step 4) Run RepEnrich2
###After we get the unique and multimap fastq, we run RepEnrich2 on the data
###Note : I have modified the original Repenrich2 script for compatibility with BioHPC
python2.7 /work/GCRB/s157114/softwares/RepEnrich2/RepEnrich2_modified.py /work/GCRB/s157114/softwares/repeatmasker_seq/mm10_repeatmasker_without_simple_repeats.txt ./ESC_H33_KO1_Input ESC_H33_KO1_Input /work/GCRB/s157114/softwares/repeatmasker_seq/setup_folder_mm10_no_simple_repeats /project/GCRB/Banaszynski_lab/shared/Yu-Ching/Repeat_Analysis/ChIP_data/LB_Nature/RepEnrich_analysis/ESC_H33_KO1_Input_R1_001_multimap.fastq /project/GCRB/Banaszynski_lab/shared/Yu-Ching/Repeat_Analysis/ChIP_data/LB_Nature/RepEnrich_analysis/ESC_H33_KO1_Input_R1_001_unique.bam --pairedend FALSE --allcountmethod TRUE

##Step 4) Run RepEnrich2
###After we get the unique and multimap fastq, we run RepEnrich2 on the data
###Note : I have modified the original Repenrich2 script for compatibility with BioHPC
python2.7 /work/GCRB/s157114/softwares/RepEnrich2/RepEnrich2_modified.py /work/GCRB/s157114/softwares/repeatmasker_seq/mm10_repeatmasker_without_simple_repeats.txt ./ESC_H33_KO2_Input ESC_H33_KO2_Input /work/GCRB/s157114/softwares/repeatmasker_seq/setup_folder_mm10_no_simple_repeats /project/GCRB/Banaszynski_lab/shared/Yu-Ching/Repeat_Analysis/ChIP_data/LB_Nature/RepEnrich_analysis/ESC_H33_KO2_Input_R1_001_multimap.fastq /project/GCRB/Banaszynski_lab/shared/Yu-Ching/Repeat_Analysis/ChIP_data/LB_Nature/RepEnrich_analysis/ESC_H33_KO2_Input_R1_001_unique.bam --pairedend FALSE --allcountmethod TRUE
