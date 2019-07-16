##Takes BAM files as input and converts to granges object.
##Plots metagenes, calculates pausing index and returns dataframe of pausing indeces for different input genes.
##Plots cumulative density pausing index graphs for WT vs KO.
##Edited by Aishwarya Gogate on 21 November, 2016.

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
legend <- c("WT", "KO")

##Add metagene function & plotting here!

#Pausing Index function
gene_list <- kgConsensus
WT_p <- pausingIndex(gene_list, WT, up = 300, down=300)
KO_p <- pausingIndex(gene_list, KO, up = 300, down=300)

write.table(WT_p, file="WT_pausing_counts.xls", sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)

#Calculate Pausing index (pause/genebody)
WT_p$PI <- WT_p$Pause/WT_p$Body
KO_p$PI <- KO_p$Pause/KO_p$Body

pause_compare <- data.frame(gene_list$gene_id,WT_p$PI,KO_p$PI)
write.table(pause_compare, file="pausing_comparison_categories_WTvsKO.xls", sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)

pause_data <- data.frame(WT_p$GeneID,WT_p$Pause,KO_p$Pause)
names(pause_data) <- c("GeneID", "WT", "KO")
pdf('boxplot_Pausing_index.pdf')
boxplot(pause_data[,2:3],outline=FALSE,border=c("black","red"), ylab="Pausing Index", names=c("WT", "KO"))
dev.off()
