#####
#Author: Rohit Setlem#

## Installing Libraries ##
install.packages("ggpubr")
library(ggpubr)
library(gplots)
##
#"Linux script: sed '1p;/yes/!d' gene_exp.diff > test.txt"
##

### Differential expression analysis ####

#Set working directory#
setwd("/Users/rohitsetlem/Desktop/Mala Lab/Mariano/1.RNA_Seq/NPvsD15/")
#Import data
data <- read.table(file.choose(),sep = "\t",header = T)
View(data)
##MA plot main
ma <- "NPvsD15"
#pdf generation
pd <- "NPvsD15_MA.pdf"
##Boxplot
bd <- "NPvsD15_Boxplot.pdf"
##Heatmap
fpkm_hm <- "NPvsD15_FPKM_hm.pdf"
fc_hm <- "NPvsD15_FC_hm.pdf"


FPKM_Analysis <- function(mydata){
  #Removing FPKM less than 1
  mydata_FPKM_1 <- mydata[which(mydata[8]>1),]
  mydata_FPKM <- mydata_FPKM_1[which(mydata_FPKM_1[9]>1),]
  #Removing any duplicate genes
  mydata_nodup <- mydata_FPKM[!duplicated(mydata_FPKM[3]),]
  #Dataframe order in ascending order
  mydata_final <- mydata_nodup[order(mydata_nodup[8]),] 
  #change 1 FC2 to 0.585 for FC 1.5 values
  Up_regulated <- mydata_final[mydata_final[10] >= 1, ] #FC 2 Up regulated genes
  Down_regulated <- mydata_final[mydata_final[10] <= -1, ] #FC 2 Down regulated genes
  #Significant genes 
  Significant <- rbind(Up_regulated,Down_regulated)
  #Saving files 
  write.csv(Up_regulated,file = paste0(ma,"_Up_regulated",".csv"),row.names = F)
  write.csv(Down_regulated,file = paste0(ma,"_Down_regulated",".csv"),row.names = F)
  write.csv(Significant,file = paste0(ma,"_Significant",".csv"),row.names = F)
  
  ## MA Plot ##
  #Taking mydata_final data 
  #selecting columns for MA plot
  MA_data <- mydata_final[,c(3,8,9,10,13)]
  rownames(MA_data) <- MA_data$gene
  #Dropping gene column
  MA_data_1 <- MA_data[2:5]
  MA_data_FPKM_mean<- MA_data_1[1:2] #selecting sample FPKM for mean
  FPKM_mean<-rowMeans(MA_data_FPKM_mean)
  #combining data
  MA_data_final <- cbind(FPKM_mean,MA_data_1[3:4])
  names<-row.names(MA_data_final)
  #change columnnames
  names(MA_data_final)[1]<-paste("baseMean")
  names(MA_data_final)[2]<-paste("log2FoldChange")
  names(MA_data_final)[3]<-paste("padj")
  MA_data_final$log2FoldChange<-as.numeric(as.character(MA_data_final$log2FoldChange))
  pdf(pd)
  maplot <- ggmaplot(MA_data_final, main = ma,
           fdr = 0.05, fc = 2, size = 0.5,legend = "top",
           genenames = as.vector(row.names(MA_data_final)),
           top=20,select.top.method = c( "fc"),
           font.label = c("bold", 11),label.rectangle = F,
           font.legend = "bold",
           font.main = "bold", 
           ggtheme = theme_classic())
  print(maplot)
  dev.off()
  
  ## Boxplot ##
  Box_data <- Significant[8:9]
  Box_data[,1:2] <- log(Box_data[1:2],2)
  pdf(bd)
  par(font.axis = 2)
  boxplot(Box_data$value_1,Box_data$value_2,main = "Boxplots",
          las = 1, names = c("Sample1","Sample2"),at =1:2,font.lab=2,whisklty = 1, lwd = 4,
          col = c("grey","darkorange"),border = "black",ylim =c(0,max(Box_data)), xlab = expression(bold("Samples")),ylab = expression(bold("Log2 FPKM")),outline=FALSE,notch = T)
  box(lwd=4)
  dev.off()
  
  ## Heatmap ##
  hm_data <- Significant[,c(3,8,9,10)]
  hm_data_log <- log(hm_data[2:3],2)
  #FPKM Heatmap
  hm_FPKM <- cbind(hm_data[1],hm_data_log)
  rownames(hm_FPKM) <- hm_FPKM$gene
  hm_FPKM_df <- hm_FPKM[2:3]
  hm_FPKM_df <- hm_FPKM_df[order(-hm_FPKM_df[1]),]
  hmcols<-colorRampPalette(c("yellow","black","blue"))
  pdf(fpkm_hm)
  heatmap.2(as.matrix(hm_FPKM_df),col = hmcols,
            symm = F,symkey=F,symbreaks=T,Rowv = FALSE,
            dendrogram='none',Colv=FALSE,breaks = c(0,1,2,3,4,5,6,7,8,9,10,11),
            margin=c(6, 6),srtCol = 45,adjCol = c(1,1),cexCol=1.2,
            xlab="Samples", ylab= "Genes",
            main="FPKM Heatmap",
            trace = "none")
  dev.off()
  #FC Heatmap
  hm_FC <- hm_data[,c(1,4)]
  hm_FC$Control <- NA
  rownames(hm_FC) <- hm_FC$gene
  hm_FC_final <- hm_FC[,c(3,2)]
  hm_FC_final_df <- hm_FC_final[order(-hm_FC_final[2]),]
  pdf(fc_hm)
  heatmap.2(as.matrix(hm_FC_final_df),col = hmcols,
            symm = F,symkey=F,symbreaks=T,Rowv = F,breaks = c(-3,-2.5,-2,-1.5,-1,-0.5,0,0.5,1,1.5,2,2.5,3),
            dendrogram='none',Colv=FALSE,na.color = "black",
            margin=c(6, 6),srtCol = 45,adjCol = c(1,1),cexCol=1.2,
            xlab="Samples", ylab= "Genes",
            main="FC Heatmap",
            trace = "none")
  dev.off()
  return(mydata)
}

FPKM_Analysis(data)
