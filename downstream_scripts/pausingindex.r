##Takes BAM files as input and converts to granges object.
##Plots metagenes, calculates pausing index and returns dataframe of pausing indeces for WT vs KO/KI for different genes.
##Plots cumulative density curves for WT vs KO/KI.
##Edited by Aishwarya Gogate on 8 November, 2016.

#load required packages
require("groHMM")
library("IRanges")
library("GenomicRanges")
library("Rsamtools")
library("GenomicAlignments")
library("GenomicFeatures")
options(mc.cores=getCores(4))

## Read alignment files
WT_rep1 <- readGAlignments("/project/GCRB/Xiaoying_lab/shared/Melodi_Analysis/downsample_bam/WT_rep1_downsampled_BAM.bam")
WT_rep2 <- readGAlignments("/project/GCRB/Xiaoying_lab/shared/Melodi_Analysis/downsample_bam/WT_rep2_downsampled_BAM.bam")
WT_rep3 <- readGAlignments("/project/GCRB/Xiaoying_lab/shared/Melodi_Analysis/downsample_bam/WT_rep3_downsampled_BAM.bam")
KI_rep1 <- readGAlignments("/project/GCRB/Xiaoying_lab/shared/Melodi_Analysis/downsample_bam/KI_rep1_downsampled_BAM.bam")
KI_rep2 <- readGAlignments("/project/GCRB/Xiaoying_lab/shared/Melodi_Analysis/downsample_bam/KI_rep2_downsampled_BAM.bam")

##convert to granges object
WT_rep1.gr <- granges(WT_rep1)
WT_rep2.gr <- granges(WT_rep2)
WT_rep3.gr <- granges(WT_rep3)
KI_rep1.gr <- granges(KI_rep1)
KI_rep2.gr <- granges(KI_rep2)

##combine replicates
WT <- c(WT_rep1.gr,WT_rep2.gr,WT_rep3.gr)
KI <- c(KI_rep1.gr, KI_rep2.gr)

##Transcripts
rf_gtf.file <- file.path("/project/GCRB/Xiaoying_lab/shared/Melodi_Analysis/downsample_bam", "ucsc_refseq_genes_mm10.gtf")
rf_txdb <- makeTxDbFromGFF(rf_gtf.file, format="gtf")
rf_kgtx <- transcripts(rf_txdb, columns=c("gene_id", "tx_id", "tx_name"))

##Collapse overlapping annotations
kgConsensus <- makeConsensusAnnotations(rf_kgtx, keytype=c("gene_id"))

##Metagene analysis (refer groHMM tutorial)
upGenes <- rf_kgtx
expReads <- mean(c(NROW(WT), NROW(KI)))
seqlevels(kgConsensus, force=TRUE) <- seqlevelsInUse(WT)
legend <- c("WT", "KI")

##Add metagene function & plotting here! (if needed)

#Pausing Index function
gene_list <- kgConsensus
WT_p <- pausingIndex(gene_list, WT, up = 300, down=300)
KI_p <- pausingIndex(gene_list, KI, up = 300, down=300)

#Calculate Pausing index (pause/genebody)
WT_p$PI <- WT_p$Pause/WT_p$Body
KI_p$PI <- KI_p$Pause/KI_p$Body

#Generate dataframe with pausing index of WT and KO/KI.
pause_index <- data.frame(WT_p$GeneID,WT_p$PI,KI_p$PI)

#Write dataframe generated in above step into an excel spreadsheet.
write.table(pause_index, file="PI_WTvsKI.xls", sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)
save.image("pausing_index_WTvsKI.RData")

##End of Script##
