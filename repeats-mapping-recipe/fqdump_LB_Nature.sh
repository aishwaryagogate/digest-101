#!/bin/tcsh

#SBATCH --job-name=fq-dump                       # job name
#SBATCH --partition=super                                 # select partion from 128GB, 256GB, 384GB, GPU and super
#SBATCH --nodes=1                                         # number of nodes requested by user
#SBATCH --time=0-96:00:00                                 # run time, format: D-H:M:S (max wallclock time)
#SBATCH --output=serialJob.%j.out                         # standard output file name
#SBATCH --error=serialJob.%j.time                         # standard error output file name
#SBATCH --mail-user=aishwarya.gogate@utsouthwestern.edu          # specify an email address
#SBATCH --mail-type=ALL                                   # send email when job status change (start, end, abortion and etc.)


##################################
#download and extract the file
###################################

#set INDEX=/project/apps_database/iGenomes/Mus_musculus/UCSC/mm10/Sequence/BWAIndex/genome.fa
set FASTQ_DUMP=/home2/s171755/sratoolkit.2.8.2-1-centos_linux64/bin/fastq-dump
set NUM_THREADS=$SLURM_CPUS_ON_NODE
set MISMATCH_PENALTY=8

module load BWA/0.7.5
module load samtools
module load picard/1.117
module load bedtools
module load sra_toolkit/2.5.1
module load UCSC_userApps
module load macs/1.4.2
module load igvtools/2.3.71


extract fastq from sra file

mkdir fastq_files

foreach sample (\
                ESC_H33WT_H3K9me3\
                ESC_H33KO1_H3K9me3\
                ESC_H33KO2_H3K9me3\
                ESC_H33WT_Input\
                ESC_H33KO1_Input\
                ESC_H33KO2_Input\
                ESC_DAXXWT_H3K9me3\
                ESC_DAXXnull_H3K9me3\
                ESC_ATRXWT_H3K9me3\
                ESC_ATRXnull_H3K9me3\
                ESC_DAXXWT_Input\
                ESC_DAXXnull_Input\
                ESC_ATRXWT_Input\
                ESC_ATRXnull_Input\
                )
  echo $sample
     extract fastq from sra file
    fastq-dump --split-3 --gzip $sample.sra

    mv *fastq.gz ./fastq_files
end
