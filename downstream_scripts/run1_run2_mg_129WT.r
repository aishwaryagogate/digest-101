#load required packages
require("groHMM")
library("IRanges")
library("GenomicRanges")
library("Rsamtools")
library("GenomicAlignments")
options(mc.cores=getCores(4))

#loading data (bam files)
WT_rep1 <- readGAlignments("/project/GCRB/Banaszynski_lab/shared/Sara/H3.3_project_Sara/GRO-Seq/run1_run2_2015_11_23_2016_02_01/sorted_bam_files/129b_13_S20.sorted.bam")
WT_rep2 <- readGAlignments("/project/GCRB/Banaszynski_lab/shared/Sara/H3.3_project_Sara/GRO-Seq/run1_run2_2015_11_23_2016_02_01/sorted_bam_files/129b_13_S11.sorted.bam")

##convert to granges object
WT_rep1.gr <- granges(WT_rep1)
WT_rep2.gr <- granges(WT_rep2)


library(GenomicFeatures)
rf_gtf.file <- file.path("/project/GCRB/Banaszynski_lab/shared/Sara/H3.3_project_Sara/GRO-Seq/run1_run2_2015_11_23_2016_02_01", "ucsc_refseq_genes_mm10.gtf")
rf_txdb <- makeTxDbFromGFF(rf_gtf.file, format="gtf")
rf_kgtx <- transcripts(rf_txdb, columns=c("gene_id", "tx_id", "tx_name"))

##Collapse overlapping annotations
kgConsensus <- makeConsensusAnnotations(rf_kgtx, keytype=c("gene_id"))

##Metagene analysis (refer groHMM tutorial)
upGenes <- rf_kgtx
expReads <- mean(c(NROW(WT_rep1.gr), NROW(WT_rep2.gr)))
seqlevels(kgConsensus, force=TRUE) <- seqlevelsInUse(WT_rep1.gr)

##Metagene around TSS
mgWT_rep1 <- runMetaGene(features=kgConsensus, reads=WT_rep1.gr, size=100, normCounts=expReads/NROW(WT_rep1.gr), sampling=FALSE, mc.cores=getOption("mc.cores"))
mgWT_rep2 <- runMetaGene(features=kgConsensus, reads=WT_rep2.gr, size=100, normCounts=expReads/NROW(WT_rep2.gr), sampling=FALSE, mc.cores=getOption("mc.cores"))

##plotting metagene
plotMetaGene <- function(POS=c(-10000:+9999), mg, MIN, MAX){
  plot(POS, mg$sense, col="red", type="h", xlim=c(-500, 13000), ylim=c(-1000,2000), ylab="Read Density", xlab="Position (relative to TSS)")
    points(POS, (-1*rev(mg$antisense)), col="blue", type="h")
    abline(mean(mg$sense[5000:8000]), 0, lty="dotted")
  }

  pdf('WT_129_rep1_genebody_metagene.pdf')
  MAX <- max(c(mgWT_rep1$sense, mgWT_rep2$sense))
  MIN <- -1*max(c(mgWT_rep1$antisense, mgWT_rep2$antisense))
  plotMetaGene(mg=mgWT_rep1, MIN=MIN, MAX=MAX)
  dev.off()

  pdf('WT_129_rep2_genebody_metagene.pdf')
  MAX <- max(c(mgWT_rep1$sense, mgWT_rep2$sense))
  MIN <- -1*max(c(mgWT_rep1$antisense, mgWT_rep2$antisense))
  plotMetaGene(mg=mgWT_rep2, MIN=MIN, MAX=MAX)
  dev.off()
