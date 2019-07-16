##Takes BAM files as input and converts to granges object.
##Plots metagenes, calculates pausing index and returns dataframe of pausing indeces for different input genes.
##Plots cumulative density pausing index graphs for WT vs KO.
##Edited by Aishwarya on 21 November, 2016.

#load required packages
require("groHMM")
library("IRanges")
library("GenomicRanges")
library("Rsamtools")
library("GenomicAlignments")
library("GenomicFeatures")
options(mc.cores=getCores(4))

## Read alignment files
WT_rep1 <- readGAlignments("/project/GCRB/Banaszynski_lab/shared/Sara/H3.3_project_Sara/GRO-Seq/run1_run2_2015_11_23_2016_02_01/sorted_bam_files/129b_13_S20.sorted.bam")
WT_rep2 <- readGAlignments("/project/GCRB/Banaszynski_lab/shared/Sara/H3.3_project_Sara/GRO-Seq/run1_run2_2015_11_23_2016_02_01/sorted_bam_files/129b_13_S11.sorted.bam")
KO_rep1 <- readGAlignments("/project/GCRB/Banaszynski_lab/shared/Sara/H3.3_project_Sara/GRO-Seq/run1_run2_2015_11_23_2016_02_01/sorted_bam_files/282_S21.sorted.bam")
KO_rep2 <- readGAlignments("/project/GCRB/Banaszynski_lab/shared/Sara/H3.3_project_Sara/GRO-Seq/run1_run2_2015_11_23_2016_02_01/sorted_bam_files/282_S12.sorted.bam")

##convert to granges object
WT_rep1.gr <- granges(WT_rep1)
WT_rep2.gr <- granges(WT_rep2)
KO_rep1.gr <- granges(KO_rep1)
KO_rep2.gr <- granges(KO_rep2)

##combine replicates
WT <- c(WT_rep1.gr,WT_rep2.gr)
KO <- c(KO_rep1.gr,KO_rep2.gr)

##Transcripts
rf_gtf.file <- file.path("/project/GCRB/Banaszynski_lab/shared/Sara/H3.3_project_Sara/GRO-Seq/run1_run2_2015_11_23_2016_02_01", "ucsc_refseq_genes_mm10.gtf")
rf_txdb <- makeTxDbFromGFF(rf_gtf.file, format="gtf")
rf_kgtx <- transcripts(rf_txdb, columns=c("gene_id", "tx_id", "tx_name"))

##Collapse overlapping annotations
kgConsensus <- makeConsensusAnnotations(rf_kgtx, keytype=c("gene_id"))

##Metagene analysis (refer groHMM tutorial)
upGenes <- rf_kgtx
expReads <- mean(c(NROW(WT), NROW(KO)))
seqlevels(kgConsensus, force=TRUE) <- seqlevelsInUse(WT)
legend <- c("WT129", "KO282")

##Add metagene function & plotting here!

#Pausing Index function
gene_list <- kgConsensus
WT_p <- pausingIndex(gene_list, WT, up = 300, down=300)
KO_p <- pausingIndex(gene_list, KO, up = 300, down=300)

#Calculate Pausing index (pause/genebody)
WT_p$PI <- WT_p$Pause/WT_p$Body
KO_p$PI <- KO_p$Pause/KO_p$Body

pause_data <- data.frame(WT_p$GeneID,WT_p$Pause,KO_p$Pause)
pause_data1 <- data.frame(WT_p$GeneID,WT_p$Body,KO_p$Body)
write.table(pause_data, file="WT_KO_pause_counts.xls", sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)
write.table(pause_data1, file="WT_KO_Body_counts.xls", sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)

#names(pause_data) <- c("GeneID", "WT", "KO")
pdf('boxplot_Pause_counts.pdf')
boxplot(pause_data[,2:3],outline=FALSE,border=c("blue","red"), ylab="Pause Counts", names=c("WT", "KO"))

pdf('boxplot_Body_counts.pdf')
boxplot(pause_data1[,2:3],outline=FALSE,border=c("blue","red"), ylab="Body Counts", names=c("WT", "KO"))

dev.off()

#ECDF (cumulative density function)
#WT_p_pause <- subset(WT_p$Pause, WT_p$Pause < max(WT_p$Pause) & WT_p$Pause > 0)
#KO_p_pause <- subset(KO_p$Pause, KO_p$Pause < max(KO_p$Pause) & KO_p$Pause > 0)
#ecdf_WT <- ecdf(log(WT_p_pause,10))
#ecdf_KO<- ecdf(log(KO_p_pause,10))
#pdf('Cumulative_density_pausing_index_WT129_KO282.pdf')
#plot(ecdf_WT, verticals=TRUE, do.points=FALSE,xlim=c(-3,4), ylab="Cumulative Fraction of Genes", xlab="Log10 Pausing Index", main="Cumulative density")
#plot(ecdf_KO, verticals=TRUE, do.points=FALSE, add=TRUE, col='red')
#legend("right", as.character(c("WT", "KO")), col=1:2,lty=1)
#dev.off()
