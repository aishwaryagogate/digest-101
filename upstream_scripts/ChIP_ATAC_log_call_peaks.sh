#!/bin/tcsh

#SBATCH --job-name=ChIP_ATAC_log                            # job name
#SBATCH --partition=super                                # select partion from 128GB, 256GB, 384GB, GPU and super
#SBATCH --nodes=1                                        # number of nodes requested by user
#SBATCH --time=0-96:00:00                                 # run time, format: D-H:M:S (max wallclock time)
#SBATCH --output=serialJob.%j.out                        # standard output file name
#SBATCH --error=serialJob.%j.time                        # standard error output file name
#SBATCH --mail-user=aishwarya.e@utsouthwestern.edu      # specify an email address
#SBATCH --mail-type=ALL                                  # send email when job status change (start, end, abortion and etc.)


#set INDEX=/project/apps_database/iGenomes/Mus_musculus/UCSC/mm10/Sequence/BWAIndex/genome.fa
#set INDEX=/project/apps_database/iGenomes/Drosophila_melanogaster/UCSC/dm3/Sequence/BWAIndex/genome.fa
set INDEX=/project/apps_database/iGenomes/Homo_sapiens/UCSC/hg19/Sequence/BWAIndex/genome.fa
set NUM_THREADS=30
set MISMATCH_PENALTY=8
#set FILE=/project/GCRB/Banaszynski_lab/shared/18_run/171009_NB501597_0173_AHMV3HBGX3/XChIP

module load BWA/0.7.5
module load samtools
module load picard/1.117
module load bedtools
module load UCSC_userApps
module load macs/1.4.2
module load igvtools/2.3.71

mkdir peaks_calls

macs14 -t ./nodup_files/filename.nodup.bam -g 1.87e9 -n ./peaks_calls/filename_peaks
macs14 -t ./nodup_files/filename.nodup.bam -g 1.87e9 -n ./peaks_calls/filename_peaks
macs14 -t ./nodup_files/filename.nodup.bam -g 1.87e9 -n ./peaks_calls/filename_peaks
macs14 -t ./nodup_files/filename.nodup.bam -g 1.87e9 -n ./peaks_calls/filename_peaks
