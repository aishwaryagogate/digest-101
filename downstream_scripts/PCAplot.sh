#!/bin/bash
#SBATCH --job-name=plotPCA                              # job name
#SBATCH --partition=super                                 # select partion from 128GB, 256GB, 384GB, GPU and super
#SBATCH --nodes=1                                         # number of nodes requested by user
#SBATCH --time=0-24:00:00                                 # run time, format: D-H:M:S (max wallclock time)
#SBATCH --output=job_PCA.%j.out                         # standard output file name
#SBATCH --error=job_PCA.%j.err                         # standard error output file name
#SBATCH --mail-user=aishwarya.gogate@utsouthwestern.edu           # specify an email address
#SBATCH --mail-type=ALL                                   # send email when job status change (start, end, abortion and etc.)

module load deeptools
##PCA plot for Addback RNA-Seq data

#multiBigwigSummary BED-file -b 282_H32_D4_RNASeq_rep1_S6_R1_001.positive.bw 282_H32_D4_RNASeq_rep2_S7_R1_001.positive.bw 282_H33_D4_RNASeq_rep1_S8_R1_001.positive.bw 282_H33_D4_RNASeq_rep2_S9_R1_001.positive.bw 282_S31A_D4_RNASeq_rep1_S10_R1_001.positive.bw 282_S31A_D4_RNASeq_rep2_S11_R1_001.positive.bw 282_S31E_D4_RNASeq_rep1_S12_R1_001.positive.bw 282_S31E_D4_RNASeq_rep2_S13_R1_001.positive.bw --BED WT129_D0_ExpressedGenes.bed --labels H3.2_1 H3.2_2 H3.3_1 H3.3_2 S31A_1 S31A_2 S31E_1 S31E_2 -out cormat_allrefseq.npz --outRawCounts rawct_allrefseq.tab
#plotPCA --corData cormat_allrefseq.npz --plotFile pca_expressed_refseq.pdf --labels H3.2_1 H3.2_2 H3.3_1 H3.3_2 S31A_1 S31A_2 S31E_1 S31E_2 --colors red red blue blue purple purple green green --plotTitle 'PCA of Addback RNAseq under expressed refseq'
