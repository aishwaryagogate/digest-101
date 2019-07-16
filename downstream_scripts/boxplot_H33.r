#load required packages
require("groHMM")
library("IRanges")
library("GenomicRanges")
library("Rsamtools")
library("GenomicAlignments")
library("GenomicFeatures")
options(mc.cores=getCores(4))
#install.packages("vioplot")
library(vioplot)

######H3K27ac-SI-rep3
WT_rep1 <- readGAlignments("/project/GCRB/Banaszynski_lab/shared/Sara/H3.3_project_Sara/ChIP-Seq/Analyis/New_analysis/nodup_normalized/129_H3K27ac_SI_rep3_norm.bam")
KO_rep1 <- readGAlignments("/project/GCRB/Banaszynski_lab/shared/Sara/H3.3_project_Sara/ChIP-Seq/Analyis/New_analysis/nodup_normalized/282_H3K27ac_SI_rep3_norm.bam")

######H3K27ac-rep1&rep2-merged
#WT_rep1 <- readGAlignments("/project/GCRB/Banaszynski_lab/shared/Sara/H3.3_project_Sara/ChIP-Seq/Analyis/New_analysis/nodup_norm_merged/129_H3K27ac_norm_merged.bam")
#KO_rep1 <- readGAlignments("/project/GCRB/Banaszynski_lab/shared/Sara/H3.3_project_Sara/ChIP-Seq/Analyis/New_analysis/nodup_norm_merged/282_H3K27ac_norm_merged.bam")

######H3K27ac-mab-rep1&rep2
#WT_rep1 <- readGAlignments("/project/GCRB/Banaszynski_lab/shared/Sara/H3.3_project_Sara/ChIP-Seq/Analyis/New_analysis/nodup_normalized/129_H3K27ac_mab_rep1_norm.bam")
#KO_rep1 <- readGAlignments("/project/GCRB/Banaszynski_lab/shared/Sara/H3.3_project_Sara/ChIP-Seq/Analyis/New_analysis/nodup_normalized/282_H3K27ac_mab_rep1_norm.bam")
#WT_rep2 <- readGAlignments("/project/GCRB/Banaszynski_lab/shared/Sara/H3.3_project_Sara/ChIP-Seq/Analyis/New_analysis/nodup_normalized/129_H3K27ac_mab_rep2_norm.bam")
#KO_rep2 <- readGAlignments("/project/GCRB/Banaszynski_lab/shared/Sara/H3.3_project_Sara/ChIP-Seq/Analyis/New_analysis/nodup_normalized/282_H3K27ac_mab_rep2_norm.bam")

######H3K64ac-rep1
#WT_rep1 <- readGAlignments("/project/GCRB/Banaszynski_lab/shared/Sara/H3.3_project_Sara/ChIP-Seq/Analyis/New_analysis/nodup_normalized/129_K64ac_rep1_norm.bam")
#KO_rep1 <- readGAlignments("/project/GCRB/Banaszynski_lab/shared/Sara/H3.3_project_Sara/ChIP-Seq/Analyis/New_analysis/nodup_normalized/282_K64ac_rep1_norm.bam")

######H3K122ac-rep1
#WT_rep1 <- readGAlignments("/project/GCRB/Banaszynski_lab/shared/Sara/H3.3_project_Sara/ChIP-Seq/Analyis/New_analysis/nodup_normalized/129_K122ac_rep1_norm.bam")
#KO_rep1 <- readGAlignments("/project/GCRB/Banaszynski_lab/shared/Sara/H3.3_project_Sara/ChIP-Seq/Analyis/New_analysis/nodup_normalized/282_K122ac_rep1_norm.bam")

######H3K9ac-rep1
#WT_rep1 <- readGAlignments("/project/GCRB/Banaszynski_lab/shared/Sara/H3.3_project_Sara/ChIP-Seq/Analyis/New_analysis/nodup_normalized/129_H3K9ac_rep1_norm.bam")
#KO_rep1 <- readGAlignments("/project/GCRB/Banaszynski_lab/shared/Sara/H3.3_project_Sara/ChIP-Seq/Analyis/New_analysis/nodup_normalized/282_H3K9ac_rep1_norm.bam")

######H3K4me1-rep1&rep2-merged
#WT_rep1 <- readGAlignments("/project/GCRB/Banaszynski_lab/shared/Sara/H3.3_project_Sara/ChIP-Seq/Analyis/New_analysis/nodup_norm_merged/129_H3K4Me1_norm_merged.bam")
#KO_rep1 <- readGAlignments("/project/GCRB/Banaszynski_lab/shared/Sara/H3.3_project_Sara/ChIP-Seq/Analyis/New_analysis/nodup_norm_merged/282_H3K4Me1_norm_merged.bam")

######H3K4me3-rep1
#WT_rep1 <- readGAlignments("/project/GCRB/Banaszynski_lab/shared/Sara/H3.3_project_Sara/ChIP-Seq/Analyis/New_analysis/nodup_normalized/129_H3K4me3_rep1_norm.bam")
#KO_rep1 <- readGAlignments("/project/GCRB/Banaszynski_lab/shared/Sara/H3.3_project_Sara/ChIP-Seq/Analyis/New_analysis/nodup_normalized/282_H3K4me3_rep1_norm.bam")

######p300-rep1
#WT_rep1 <- readGAlignments("/project/GCRB/Banaszynski_lab/shared/Sara/H3.3_project_Sara/ChIP-Seq/Analyis/New_analysis/nodup_normalized/129_P300_III_rep1_norm.bam")
#KO_rep1 <- readGAlignments("/project/GCRB/Banaszynski_lab/shared/Sara/H3.3_project_Sara/ChIP-Seq/Analyis/New_analysis/nodup_normalized/282_P300_III_rep1_norm.bam")

######p300-rep2
#WT_rep1 <- readGAlignments("/project/GCRB/Banaszynski_lab/shared/Sara/H3.3_project_Sara/ChIP-Seq/Analyis/New_analysis/nodup_normalized/129_p300_IV_rep2_norm.bam")
#KO_rep1 <- readGAlignments("/project/GCRB/Banaszynski_lab/shared/Sara/H3.3_project_Sara/ChIP-Seq/Analyis/New_analysis/nodup_normalized/282_p300_IV_rep2_norm.bam")

######ATAC-old-rep1
#WT_rep1 <- readGAlignments("/project/GCRB/Banaszynski_lab/shared/Sara/H3.3_project_Sara/ChIP-Seq/Analyis/New_analysis/nodup_normalized/old129b_rep1_norm.bam")
#KO_rep1 <- readGAlignments("/project/GCRB/Banaszynski_lab/shared/Sara/H3.3_project_Sara/ChIP-Seq/Analyis/New_analysis/nodup_normalized/old282_rep1_norm.bam")

######ATAC-new-rep2
#WT_rep1 <- readGAlignments("/project/GCRB/Banaszynski_lab/shared/Sara/H3.3_project_Sara/ChIP-Seq/Analyis/New_analysis/nodup_normalized/new129b_rep2_norm.bam")
#KO_rep1 <- readGAlignments("/project/GCRB/Banaszynski_lab/shared/Sara/H3.3_project_Sara/ChIP-Seq/Analyis/New_analysis/nodup_normalized/new282_rep2_norm.bam")

######p300Ac-rep2
#WT_rep1 <- readGAlignments("/project/GCRB/Banaszynski_lab/shared/Sara/H3.3_project_Sara/ChIP-Seq/Analyis/New_analysis/nodup_normalized/129_p300Ac_rep2_norm.bam")
#KO_rep1 <- readGAlignments("/project/GCRB/Banaszynski_lab/shared/Sara/H3.3_project_Sara/ChIP-Seq/Analyis/New_analysis/nodup_normalized/282_p300Ac_rep2_norm.bam")

######p300Ac-rep3
#WT_rep1 <- readGAlignments("/project/GCRB/Banaszynski_lab/shared/Sara/H3.3_project_Sara/ChIP-Seq/Analyis/New_analysis/nodup_normalized/129_p300Ac_rep3_norm.bam")
#KO_rep1 <- readGAlignments("/project/GCRB/Banaszynski_lab/shared/Sara/H3.3_project_Sara/ChIP-Seq/Analyis/New_analysis/nodup_normalized/282_p300Ac_rep3_norm.bam")

######p300-rep3
#WT_rep1 <- readGAlignments("/project/GCRB/Banaszynski_lab/shared/Sara/H3.3_project_Sara/ChIP-Seq/Analyis/New_analysis/nodup_normalized/129_p300_V_rep3_norm.bam")
#KO_rep1 <- readGAlignments("/project/GCRB/Banaszynski_lab/shared/Sara/H3.3_project_Sara/ChIP-Seq/Analyis/New_analysis/nodup_normalized/282_p300_V_rep3_norm.bam")

######ATAC-rep3
#WT_rep1 <- readGAlignments("/project/GCRB/Banaszynski_lab/shared/Sara/H3.3_project_Sara/ChIP-Seq/Analyis/New_analysis/nodup_normalized/129_III_ATAC_rep3_norm.bam")
#KO_rep1 <- readGAlignments("/project/GCRB/Banaszynski_lab/shared/Sara/H3.3_project_Sara/ChIP-Seq/Analyis/New_analysis/nodup_normalized/282_III_ATAC_rep3_norm.bam")

######p300_mAb_rep4
#WT_rep1 <- readGAlignments("/project/GCRB/Banaszynski_lab/shared/Sara/H3.3_project_Sara/ChIP-Seq/Analyis/New_analysis/nodup_normalized/129_p300_mAb_rep4_norm.bam")
#KO_rep1 <- readGAlignments("/project/GCRB/Banaszynski_lab/shared/Sara/H3.3_project_Sara/ChIP-Seq/Analyis/New_analysis/nodup_normalized/282_p300_mAb_rep4_norm.bam")

######ATAC_D4_rep1
#WT_rep1 <- readGAlignments("/project/GCRB/Banaszynski_lab/shared/Sara/H3.3_project_Sara/ChIP-Seq/Analyis/New_analysis/nodup_normalized/129_D4_ATAC_rep1_norm.bam")
#KO_rep1 <- readGAlignments("/project/GCRB/Banaszynski_lab/shared/Sara/H3.3_project_Sara/ChIP-Seq/Analyis/New_analysis/nodup_normalized/282_D4_ATAC_rep1_norm.bam")

######H3K27ac_D4_129rep3_282rep2
#WT_rep1 <- readGAlignments("/project/GCRB/Banaszynski_lab/shared/Sara/H3.3_project_Sara/ChIP-Seq/Analyis/New_analysis/nodup_normalized/129_H3K27ac_D4_rep3_norm.bam")
#KO_rep1 <- readGAlignments("/project/GCRB/Banaszynski_lab/shared/Sara/H3.3_project_Sara/ChIP-Seq/Analyis/New_analysis/nodup_normalized/282_H3K27ac_D4_rep2_norm.bam")

WT_rep1.gr <- granges(WT_rep1)
#WT_rep2.gr <- granges(WT_rep2)
KO_rep1.gr <- granges(KO_rep1)
#KO_rep2.gr <- granges(KO_rep2)

##Assign to variable
WT <- WT_rep1.gr
KO <- KO_rep1.gr

##Combine replicates
##If there are replicates, only then.
#WT <- c(WT_rep1.gr,WT_rep2.gr)
#KO <- c(KO_rep1.gr,KO_rep2.gr)

#take either enhancer or refseq.bed
#enhancers <- import("p300_narrowpeaks.bed", format = "BED")
enhancers <- import("p300merged_2500bp_distal_peaks.bed", format = "BED")
#superenhancers <- import("super_enhancers_Whyte_mm10.bed", format = "BED")
#genes <- import("mm10_RefSeqgenes.bed", format = "BED")

## Find average
library_WT_129 <- NROW(WT)
library_KO_282 <- NROW(KO)

# Calculate RPKM
#Under enhancers
## Import files RPKM =   numReads / ( geneLength/1000 * totalNumReads/1,000,000 )
rpkm_WT_129 <- countOverlaps(enhancers, WT) / (width(enhancers)/1000 * library_WT_129/1000000)
rpkm_KO_282 <- countOverlaps(enhancers, KO) / (width(enhancers)/1000 * library_KO_282/1000000)
rpkm_enhancers_p300 <- data.frame(rpkm_WT_129,rpkm_KO_282)
head(rpkm_enhancers_p300)
write.table(rpkm_enhancers_p300, file="rpkm_H3K27ac_SI_rep3_under_ALLenhancers.xls", sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)

#pdf('boxplot_p300_rep1_rep2_merged_distalDOWN.pdf')
#boxplot(rpkm_enhancers_p300, col=(c("blue","red3")), main="p300 merged under distal (DOWN)", ylab="RPKM",outline=FALSE, notch=TRUE, names=c("WT", "H3.3 KO"),yaxt="n", cex.axis=1,las=2,lwd=4,lty=1)
#axis(2, at=seq(0,12,.5))
#axis(2, at=seq(0,1.5,.2))
#dev.off()

#wilcox.test(rpkm_enhancers_p300$rpkm_WT_129, rpkm_enhancers_p300$rpkm_KO_282,conf.int=TRUE)

##End of Script##
