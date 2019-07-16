#load required packages
require("groHMM")
library("IRanges")
library("GenomicRanges")
library("Rsamtools")
library("GenomicAlignments")
library("GenomicFeatures")
options(mc.cores=getCores(4))

## Read alignment files
#ES_WT_rep2 <- readGAlignments("/project/GCRB/Xiaoying_lab/shared/Melodi_NGSdatasets/GROseq/FirstRound/Run38_3_MT_Sample3-35037/align-bwa.sh-1.0.0/Run38-3-MT-Sample3.sorted.bam")
#ES_WT_rep3 <- readGAlignments("/project/GCRB/Xiaoying_lab/shared/Melodi_NGSdatasets/GROseq/FirstRound/Run38_6_MT_Sample6-35040/align-bwa.sh-1.0.0/Run38-6-MT-Sample6.sorted.bam")
#WT_rep1 <- readGAlignments("/project/GCRB/Xiaoying_lab/shared/Melodi_Analysis/downsample_bam/WT_rep1_downsampled_BAM.bam")
#WT_rep2 <- readGAlignments("/project/GCRB/Xiaoying_lab/shared/Melodi_Analysis/downsample_bam/WT_rep2_downsampled_BAM.bam")
#WT_rep3 <- readGAlignments("/project/GCRB/Xiaoying_lab/shared/Melodi_Analysis/downsample_bam/WT_rep3_downsampled_BAM.bam")

##convert to granges object
#ES_WT_rep2.gr <- granges(ES_WT_rep2)
#ES_WT_rep3.gr <- granges(ES_WT_rep3)
#WT_rep1.gr <- granges(WT_rep1)
#WT_rep2.gr <- granges(WT_rep2)
#WT_rep3.gr <- granges(WT_rep3)

##combine replicates
#OR <- c(ES_WT_rep2.gr,ES_WT_rep3.gr)
#WT <- c(WT_rep1.gr,WT_rep2.gr,WT_rep3.gr)

## Read alignment files

WT_rep1 <- readGAlignments("/project/GCRB/Xiaoying_lab/shared/Xiuli_Analysis/downsample_bam/CD34_con1_rep1_downsampled.bam")
WT_rep2 <- readGAlignments("/project/GCRB/Xiaoying_lab/shared/Xiuli_Analysis/downsample_bam/CD34_con2_rep2_downsampled.bam")
WT_rep3 <- readGAlignments("/project/GCRB/Xiaoying_lab/shared/Xiuli_Analysis/downsample_bam/CD34_con3_rep3_downsampled.bam")
#KO_rep1 <- readGAlignments("/project/GCRB/Xiaoying_lab/shared/Xiuli_Analysis/downsample_bam/NelfB_kD_rep1_downsampled.bam")
#KO_rep2 <- readGAlignments("/project/GCRB/Xiaoying_lab/shared/Xiuli_Analysis/downsample_bam/NelfB_kD_rep2_downsampled.bam")
#KO_rep3 <- readGAlignments("/project/GCRB/Xiaoying_lab/shared/Xiuli_Analysis/downsample_bam/NelfB_kD_rep3_downsampled.bam")
KD_rep1 <- readGAlignments("/project/GCRB/Xiaoying_lab/shared/Xiuli_Analysis/downsample_bam/NelfE_kD_rep1_downsampled.bam")
KD_rep2 <- readGAlignments("/project/GCRB/Xiaoying_lab/shared/Xiuli_Analysis/downsample_bam/NelfE_kD_rep2_downsampled.bam")
KD_rep3 <- readGAlignments("/project/GCRB/Xiaoying_lab/shared/Xiuli_Analysis/downsample_bam/NelfE_kD_rep3_downsampled.bam")
#KI_rep1 <- readGAlignments("/project/GCRB/Xiaoying_lab/shared/Xiuli_Analysis/downsample_bam/MYE_d1_rep1_downsampled.bam")
#KI_rep2 <- readGAlignments("/project/GCRB/Xiaoying_lab/shared/Xiuli_Analysis/downsample_bam/MYE_d1_rep2_downsampled.bam")

##convert to granges object
WT_rep1.gr <- granges(WT_rep1)
WT_rep2.gr <- granges(WT_rep2)
WT_rep3.gr <- granges(WT_rep3)
#KO_rep1.gr <- granges(KO_rep1)
#KO_rep2.gr <- granges(KO_rep2)
#KO_rep3.gr <- granges(KO_rep3)
KD_rep1.gr <- granges(KD_rep1)
KD_rep2.gr <- granges(KD_rep2)
KD_rep3.gr <- granges(KD_rep3)
#KI_rep1.gr <- granges(KI_rep1)
#KI_rep2.gr <- granges(KI_rep2)

##combine replicates
#WT <- c(WT_rep1.gr,WT_rep2.gr,WT_rep3.gr)
#KO <- c(KO_rep1.gr,KO_rep2.gr,KO_rep3.gr)
#KD <- c(KD_rep1.gr,KD_rep2.gr,KD_rep3.gr)
#KI <- c(KI_rep1.gr, KI_rep2.gr)

##combine replicates
WT <- c(WT_rep1.gr,WT_rep2.gr,WT_rep3.gr)
OR <- c(KD_rep1.gr,KD_rep2.gr,KD_rep3.gr)


##Transcripts
rf_gtf.file <- file.path("/project/GCRB/Xiaoying_lab/shared/Xiuli_Analysis/downsample_bam", "ucsc_refseq_genes.gtf")
rf_txdb <- makeTxDbFromGFF(rf_gtf.file, format="gtf")
rf_kgtx <- transcripts(rf_txdb, columns=c("gene_id", "tx_id", "tx_name"))

##Collapse overlapping annotations
kgConsensus <- makeConsensusAnnotations(rf_kgtx, keytype=c("gene_id"))

##Metagene analysis (refer groHMM tutorial)
upGenes <- rf_kgtx
expReads <- mean(c(NROW(WT), NROW(OR)))
seqlevels(kgConsensus, force=TRUE) <- seqlevelsInUse(WT)
#legend <- c("WT", "Mye-d1")

##Add metagene function & plotting here!

#Pausing Index function
gene_list <- kgConsensus
WT_p <- pausingIndex(gene_list, WT, up = 300, down=300)
OR_p <- pausingIndex(gene_list, OR, up = 300, down=300)

#Calculate Pausing index (pause/genebody)
WT_p$PI <- WT_p$Pause/WT_p$Body
OR_p$PI <- OR_p$Pause/OR_p$Body

#ECDF (cumulative density function)
WT_p_pause <- subset(WT_p$Pause, WT_p$Pause < max(WT_p$Pause) & WT_p$Pause > 0)
OR_p_pause <- subset(OR_p$Pause, OR_p$Pause < max(OR_p$Pause) & OR_p$Pause > 0)
ecdf_WT <- ecdf(log(WT_p_pause,2))
ecdf_OR<- ecdf(log(OR_p_pause,2))
pdf('Cumulative_density_log2PI_WT_NelfEKD_CORRECTCOLORS.pdf')
plot(ecdf_WT, verticals=TRUE, do.points=FALSE,xlim=c(-6,8), col='red', ylab="Cumulative Fraction of Genes", xlab="Log2 Pausing Index", main="Cumulative density")
plot(ecdf_OR, verticals=TRUE, do.points=FALSE, add=TRUE, col='blue')
#legend("right", as.character(c("WT", "Mye-d1")), col=1:2,lty=1)
dev.off()
