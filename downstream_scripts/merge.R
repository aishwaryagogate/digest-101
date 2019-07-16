##Code to join two data files by a specific column.
##Last edit by Aishwarya on 8-July-2019.

#Read in your reference genes file
info<-read.table(file="NNG_diff_ATAC_peaks.txt", sep="\t", header=T)
head(info)

#Read in your RNA-seq DE gene counts
refinfo<-read.table(file="RNAseq_DE_genes.txt", sep="\t", header=T)
head(refinfo)

#Merge by the gene ID column
total<-merge(info, refinfo, by="ID",all.x = TRUE)
write.table(total, file="output.txt", sep="\t")

#Note: Make sure that the header of the column containing the gene names is "ID" in both your input files, as you are merging using that identifier.
