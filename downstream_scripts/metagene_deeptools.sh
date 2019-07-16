#!/bin/bash
#SBATCH --job-name=TSS-to-TES                              # job name
#SBATCH --partition=super                                 # select partion from 128GB, 256GB, 384GB, GPU and super
#SBATCH --nodes=1                                         # number of nodes requested by user
#SBATCH --time=0-24:00:00                                 # run time, format: D-H:M:S (max wallclock time)
#SBATCH --output=job_prof_HM.%j.out                         # standard output file name
#SBATCH --error=job_prof_HM.%j.err                         # standard error output file name
#SBATCH --mail-user=aishwarya.gogate@utsouthwestern.edu           # specify an email address
#SBATCH --mail-type=ALL                                   # send email when job status change (start, end, abortion and etc.)

module load deeptools

##WT & H3.3 KO H3K27ac metagene under Refseq genes (TSS to TES)
computeMatrix scale-regions --beforeRegionStartLength 1000 --regionBodyLength 5000 --afterRegionStartLength 1000 -R mm10_RefSeqgenes.bed --sortRegions no -S /project/GCRB/Banaszynski_lab/shared/Sara/H3.3_project_Sara/ChIP-Seq/Analyis/New_analysis/bigwig_normalized/129_H3K27ac_SI_rep3_norm.nodup.bw /project/GCRB/Banaszynski_lab/shared/Sara/H3.3_project_Sara/ChIP-Seq/Analyis/New_analysis/bigwig_normalized/282_H3K27ac_SI_rep3_norm.nodup.bw --skipZeros -o matrix_TSStoTES.gz --outFileSortedRegions matrix_TSStoTES.bed --numberOfProcessors max
plotHeatmap -m matrix_TSStoTES.gz -out 129_282_H3K27ac_metagene_underRefseqGenes_TSStoTES_heatmap.pdf --heatmapHeight 10 --refPointLabel center --samplesLabel 129 282 --colorMap Blues --plotTitle '129, 282 H3K27ac under Refseq' --dpi 300 --whatToShow 'heatmap and colorbar' --boxAroundHeatmaps no
plotProfile -m matrix_TSStoTES.gz -out 129_282_H3K27ac_metagene_underRefseqGenes_TSStoTES_profile.pdf --perGroup --color blue red --samplesLabel 129 282 --startLabel TSS --plotTitle "129, 282 H3K27ac under Refseq"

##Note: Can be plotted under enhancers/peaks file too.
