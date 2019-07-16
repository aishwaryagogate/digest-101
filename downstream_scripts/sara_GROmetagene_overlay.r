##Takes BAM files as input and converts to granges object.
##Generates an overlay metagene plot of WT vs KO.
##Edited by Aishwarya Gogate on 9 February, 2017.

#load required packages
require("groHMM")
library("IRanges")
library("GenomicRanges")
library("Rsamtools")
library("GenomicAlignments")
library("GenomicFeatures")
options(mc.cores=getCores(4))

## Read alignment files
WT_rep1 <- readGAlignments("/project/GCRB/Banaszynski_lab/shared/Sara/H3.3_project_Sara/GRO-Seq/run1_run2_2015_11_23_2016_02_01/sorted_bam_files/129b_13_S11.sorted.bam")
WT_rep2 <- readGAlignments("/project/GCRB/Banaszynski_lab/shared/Sara/H3.3_project_Sara/GRO-Seq/run1_run2_2015_11_23_2016_02_01/sorted_bam_files/129b_13_S20.sorted.bam")
KO_rep1 <- readGAlignments("/project/GCRB/Banaszynski_lab/shared/Sara/H3.3_project_Sara/GRO-Seq/run1_run2_2015_11_23_2016_02_01/sorted_bam_files/282_S12.sorted.bam")
KO_rep2 <- readGAlignments("/project/GCRB/Banaszynski_lab/shared/Sara/H3.3_project_Sara/GRO-Seq/run1_run2_2015_11_23_2016_02_01/sorted_bam_files/282_S21.sorted.bam")

WT_rep1.gr <- granges(WT_rep1)
WT_rep2.gr <- granges(WT_rep2)
KO_rep1.gr <- granges(KO_rep1)
KO_rep2.gr <- granges(KO_rep2)

##combine replicates
WT <- c(WT_rep1.gr,WT_rep2.gr)
KO <- c(KO_rep1.gr,KO_rep2.gr)

##Transcripts
rf_gtf.file <- file.path("/project/GCRB/Banaszynski_lab/shared/Sara/H3.3_project_Sara/GRO-Seq/run1_run2_2015_11_23_2016_02_01/rpkm_plots", "ucsc_refseq_genes_mm10.gtf")
rf_txdb <- makeTxDbFromGFF(rf_gtf.file, format="gtf")
rf_kgtx <- transcripts(rf_txdb, columns=c("gene_id", "tx_id", "tx_name"))

##Collapse overlapping annotations
kgConsensus <- makeConsensusAnnotations(rf_kgtx, keytype=c("gene_id"))

##Metagene analysis (refer groHMM tutorial)
upGenes <- rf_kgtx
#expReads <- mean(c(NROW(WT), NROW(KO)))
#Since the y-axis becomes a normalized read density between WT&KI and WT&KO, its range on the metagene in both plots is inconsistent.
#Hence using reads per million (RPM) approach where expReads= 1 million.
expReads <- 1000000
legend <- c("WT", "KO")
seqlevels(kgConsensus, force=TRUE) <- seqlevelsInUse(WT)

##plotting metagene
plotMetaGene <- function(filename, mg_a, mg_b, MIN, MAX, legend){
    pdf(filename)
    POS=c(-10000:+9999)
    print(paste(MIN, MAX))
    par(font=3, font.lab=2, font.axis=2, mar=c(5.5, 5.5, 2, 2) + 0.3)
    plot(POS, mg_a$sense, col="red", type="l", lwd=3, xlim=c(-6000, 6000), cex.axis=2, cex.lab=2, ylim=c(MIN, MAX), ylab="RPM", xlab="Position (relative to TSS)")
    lines(POS, (-1*rev(mg_a$antisense)), col="red", type="l")
    lines(POS, mg_b$sense, col="blue", type="l", lwd=3, xlim=c(-6000, 6000))
    lines(POS, (-1*rev(mg_b$antisense)), col="blue", type="l")
    lines(POS,0*POS)
    legend("topright", legend=legend, lty=1, lwd=3, col=c("red", "blue"))
    dev.off()
}

##Metagene around TSS
mg_WT <- runMetaGene(features=kgConsensus, reads=WT, size=100, normCounts=expReads/NROW(WT), sampling=FALSE)
mg_KO<- runMetaGene(features=kgConsensus, reads=KO, size=100, normCounts=expReads/NROW(KO), sampling=FALSE)
MAX <- max(c(mg_WT$sense, mg_KO$sense))
MIN <- -1*max(c(mg_WT$antisense, mg_KO$antisense))

plotMetaGene("WT_vs_KO_RPMoverlay_6kb.pdf", mg_WT, mg_KO, MIN, MAX, legend)

##End of Script.
