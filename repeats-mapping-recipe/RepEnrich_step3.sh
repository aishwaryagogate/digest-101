#!/bin/bash
#SBATCH -J RepEnrich_step3                           			 # job name
#SBATCH -p super                                 # select partion from 128GB, 256GB, 384GB, GPU and super
#SBATCH -N 1                                     			 # number of nodes requested by user
#SBATCH -t  0-96:00:00                                 # run time, format: D-H:M:S (max wallclock time)
#SBATCH --output=RepEnrich_step3.%j.out                         # standard output file name
#SBATCH --error=RepEnrich_step3.%j.err                         # standard error output file name
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

##Step 3: Map the data to the genome
##LB Nature data H3K9me3 ChIP-Seq for WT, KO1 and KO2
##Command Usage##
#gunzip -c /project/GCRB/Banaszynski_lab/shared/Yu-Ching/Repeat_Analysis/ChIP_data/LB_Nature/fastq_files/Sample_R1_001.fastq.gz > Sample_R1_001.fastq
#gunzip -c /project/GCRB/Banaszynski_lab/shared/Yu-Ching/Repeat_Analysis/ChIP_data/LB_Nature/fastq_files/Sample_R2_001.fastq.gz > Sample_R2_001.fastq
#bowtie2 -q -t -x /project/apps_database/iGenomes/Mus_musculus/UCSC/mm10/Sequence/Bowtie2Index/genome -1 Sample_R1_001.fastq -2 Sample_R2_001.fastq -S Sample.sam
#Convert SAM to BAM
#samtools view -bS Sample.sam > Sample.bam
#RUN the RepEnrich2_subset.py to output discrete files for uniquely and multi-mapping reads
#python2.7 /work/GCRB/s157114/softwares/RepEnrich2/RepEnrich2_subset.py Sample.bam 30 Sample --pairedend TRUE


##Note LB Nature data is single-end. So modify above commands.
##First trying for H3.3 WT, KO1 & KO2 H3K9me3 ChIP AND their Inputs only!!

####################WT and it's Input################
#gunzip -c /project/GCRB/Banaszynski_lab/shared/Yu-Ching/Repeat_Analysis/ChIP_data/LB_Nature/fastq_files/ESC_H33_WT_H3K9me3_R1_001.fastq.gz > ESC_H33_WT_H3K9me3_R1_001.fastq
#bowtie2 -q -t -x /project/apps_database/iGenomes/Mus_musculus/UCSC/mm10/Sequence/Bowtie2Index/genome -U ESC_H33_WT_H3K9me3_R1_001.fastq -S ESC_H33_WT_H3K9me3_R1_001.sam
#Convert SAM to BAM
#samtools view -bS ESC_H33_WT_H3K9me3_R1_001.sam > ESC_H33_WT_H3K9me3_R1_001.bam
#RUN the RepEnrich2_subset.py to output discrete files for uniquely and multi-mapping reads
#python2.7 /work/GCRB/s157114/softwares/RepEnrich2/RepEnrich2_subset.py ESC_H33_WT_H3K9me3_R1_001.bam 30 ESC_H33_WT_H3K9me3_R1_001 --pairedend FALSE

#gunzip -c /project/GCRB/Banaszynski_lab/shared/Yu-Ching/Repeat_Analysis/ChIP_data/LB_Nature/fastq_files/ESC_H33_WT_Input_R1_001.fastq.gz > ESC_H33_WT_Input_R1_001.fastq
#bowtie2 -q -t -x /project/apps_database/iGenomes/Mus_musculus/UCSC/mm10/Sequence/Bowtie2Index/genome -U ESC_H33_WT_Input_R1_001.fastq -S ESC_H33_WT_Input_R1_001.sam
#Convert SAM to BAM
#samtools view -bS ESC_H33_WT_Input_R1_001.sam > ESC_H33_WT_Input_R1_001.bam
#RUN the RepEnrich2_subset.py to output discrete files for uniquely and multi-mapping reads
#python2.7 /work/GCRB/s157114/softwares/RepEnrich2/RepEnrich2_subset.py ESC_H33_WT_Input_R1_001.bam 30 ESC_H33_WT_Input_R1_001 --pairedend FALSE


####################KO1, KO2 and their respective Inputs###############
gunzip -c /project/GCRB/Banaszynski_lab/shared/Yu-Ching/Repeat_Analysis/ChIP_data/LB_Nature/fastq_files/ESC_H33_KO1_H3K9me3_R1_001.fastq.gz > ESC_H33_KO1_H3K9me3_R1_001.fastq
bowtie2 -q -t -x /project/apps_database/iGenomes/Mus_musculus/UCSC/mm10/Sequence/Bowtie2Index/genome -U ESC_H33_KO1_H3K9me3_R1_001.fastq -S ESC_H33_KO1_H3K9me3_R1_001.sam
#Convert SAM to BAM
samtools view -bS ESC_H33_KO1_H3K9me3_R1_001.sam > ESC_H33_KO1_H3K9me3_R1_001.bam
#RUN the RepEnrich2_subset.py to output discrete files for uniquely and multi-mapping reads
python2.7 /work/GCRB/s157114/softwares/RepEnrich2/RepEnrich2_subset.py ESC_H33_KO1_H3K9me3_R1_001.bam 30 ESC_H33_KO1_H3K9me3_R1_001 --pairedend FALSE

gunzip -c /project/GCRB/Banaszynski_lab/shared/Yu-Ching/Repeat_Analysis/ChIP_data/LB_Nature/fastq_files/ESC_H33_KO2_H3K9me3_R1_001.fastq.gz > ESC_H33_KO2_H3K9me3_R1_001.fastq
bowtie2 -q -t -x /project/apps_database/iGenomes/Mus_musculus/UCSC/mm10/Sequence/Bowtie2Index/genome -U ESC_H33_KO2_H3K9me3_R1_001.fastq -S ESC_H33_KO2_H3K9me3_R1_001.sam
#Convert SAM to BAM
samtools view -bS ESC_H33_KO2_H3K9me3_R1_001.sam > ESC_H33_KO2_H3K9me3_R1_001.bam
#RUN the RepEnrich2_subset.py to output discrete files for uniquely and multi-mapping reads
python2.7 /work/GCRB/s157114/softwares/RepEnrich2/RepEnrich2_subset.py ESC_H33_KO2_H3K9me3_R1_001.bam 30 ESC_H33_KO2_H3K9me3_R1_001 --pairedend FALSE

gunzip -c /project/GCRB/Banaszynski_lab/shared/Yu-Ching/Repeat_Analysis/ChIP_data/LB_Nature/fastq_files/ESC_H33_KO1_Input_R1_001.fastq.gz > ESC_H33_KO1_Input_R1_001.fastq
bowtie2 -q -t -x /project/apps_database/iGenomes/Mus_musculus/UCSC/mm10/Sequence/Bowtie2Index/genome -U ESC_H33_KO1_Input_R1_001.fastq -S ESC_H33_KO1_Input_R1_001.sam
#Convert SAM to BAM
samtools view -bS ESC_H33_KO1_Input_R1_001.sam > ESC_H33_KO1_Input_R1_001.bam
#RUN the RepEnrich2_subset.py to output discrete files for uniquely and multi-mapping reads
python2.7 /work/GCRB/s157114/softwares/RepEnrich2/RepEnrich2_subset.py ESC_H33_KO1_Input_R1_001.bam 30 ESC_H33_KO1_Input_R1_001 --pairedend FALSE

gunzip -c /project/GCRB/Banaszynski_lab/shared/Yu-Ching/Repeat_Analysis/ChIP_data/LB_Nature/fastq_files/ESC_H33_KO2_Input_R1_001.fastq.gz > ESC_H33_KO2_Input_R1_001.fastq
bowtie2 -q -t -x /project/apps_database/iGenomes/Mus_musculus/UCSC/mm10/Sequence/Bowtie2Index/genome -U ESC_H33_KO2_Input_R1_001.fastq -S ESC_H33_KO2_Input_R1_001.sam
#Convert SAM to BAM
samtools view -bS ESC_H33_KO2_Input_R1_001.sam > ESC_H33_KO2_Input_R1_001.bam
#RUN the RepEnrich2_subset.py to output discrete files for uniquely and multi-mapping reads
python2.7 /work/GCRB/s157114/softwares/RepEnrich2/RepEnrich2_subset.py ESC_H33_KO2_Input_R1_001.bam 30 ESC_H33_KO2_Input_R1_001 --pairedend FALSE
