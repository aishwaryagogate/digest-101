#!/bin/tcsh

#SBATCH --job-name=ftp                           # job name
#SBATCH --partition=super                                 # select partion from 128GB, 256GB, 384GB, GPU and super
#SBATCH --nodes=1                                         # number of nodes requested by user
#SBATCH --time=0-96:00:00                                 # run time, format: D-H:M:S (max wallclock time)
#SBATCH --output=serialJob.%j.out                         # standard output file name
#SBATCH --error=serialJob.%j.time                         # standard error output file name
#SBATCH --mail-user=aishwarya.gogate@utsouthwestern.edu          # specify an email address
#SBATCH --mail-type=ALL                                   # send email when job status change (start, end, abortion and etc.)


##################################
#
#download and extract the file, mapping and remove duplicates
#Please use this script as of 2019
#
###################################
set INDEX=/project/apps_database/iGenomes/Mus_musculus/UCSC/mm10/Sequence/BWAIndex/genome.fa
set NUM_THREADS=$SLURM_CPUS_ON_NODE
set MISMATCH_PENALTY=8

module load sra_toolkit/2.5.1
module load BWA/0.7.5
module load samtools
module load picard/1.117
module load bedtools
module load UCSC_userApps

#download the data
##Usage: FTP_path_to_GEO_sample/file.sra > renamed_file.sra

wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR150/SRR1508933/SRR1508933.sra > ESC_H33WT_H3K9me3.sra
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR150/SRR1508934/SRR1508934.sra > ESC_H33KO1_H3K9me3.sra
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR150/SRR1508935/SRR1508935.sra > ESC_H33KO2_H3K9me3.sra
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR150/SRR1508942/SRR1508942.sra > ESC_H33WT_Input.sra
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR150/SRR1508943/SRR1508943.sra > ESC_H33KO1_Input.sra
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR150/SRR1508944/SRR1508944.sra > ESC_H33KO2_Input.sra
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR201/SRR2014810/SRR2014810.sra > ESC_DAXXWT_H3K9me3.sra
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR201/SRR2014811/SRR2014811.sra > ESC_DAXXnull_H3K9me3.sra
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR201/SRR2014819/SRR2014819.sra > ESC_ATRXWT_H3K9me3.sra
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR201/SRR2014820/SRR2014820.sra > ESC_ATRXnull_H3K9me3.sra
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR201/SRR2014812/SRR2014812.sra > ESC_DAXXWT_Input.sra
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR150/SRR1508949/SRR1508949.sra > ESC_DAXXnull_Input.sra
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR201/SRR2014823/SRR2014823.sra > ESC_ATRXWT_Input.sra
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR201/SRR2014824/SRR2014824.sra > ESC_ATRXnull_Input.sra
