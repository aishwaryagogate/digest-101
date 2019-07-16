##REPENRICH2 STEPS
##FOLLOW THE STEPS AVAILABLE HERE
https://github.com/nerettilab/RepEnrich2


##Step 1) Attain repetitive element annotation
#Download the organism specific repeatmasker file from repeatmasker.org
#http://www.repeatmasker.org/genomicDatasets/RMGenomicDatasets.html
#For instance, find the hg19.fa.out.gz/mm10.fa.out.gz and download.
#Once you have downloaded the file you can unzip it and rename it
gunzip mm10.fa.out.gz
mv mm10.fa.out mm10_repeatmasker.txt


##Step 2) Run the setup for RepEnrich2
module load iGenomes/2013-03-25
module load bedtools/2.26.0
module load bowtie2/2.2.8-intel
module load samtools/1.6
module load perl/5.28.0
python2.7 /work/GCRB/s157114/softwares/RepEnrich2/RepEnrich2_setup.py hg19_repeatmasker.txt /project/apps_database/iGenomes/Mus_musculus/UCSC/mm10/Sequence/Bowtie2Index/genome.fa ./setup_folder_mm10/


##Step 3) Map the data to the genome
##Command Usage##
#gunzip -c /project/GCRB/Banaszynski_lab/shared/Yu-Ching/Repeat_Analysis/ChIP_data/LB_Nature/fastq_files/Sample_R1_001.fastq.gz > Sample_R1_001.fastq
#gunzip -c /project/GCRB/Banaszynski_lab/shared/Yu-Ching/Repeat_Analysis/ChIP_data/LB_Nature/fastq_files/Sample_R2_001.fastq.gz > Sample_R2_001.fastq
#bowtie2 -q -t -x /project/apps_database/iGenomes/Mus_musculus/UCSC/mm10/Sequence/Bowtie2Index/genome -1 Sample_R1_001.fastq -2 Sample_R2_001.fastq -S Sample.sam
#Convert SAM to BAM
#samtools view -bS Sample.sam > Sample.bam
#RUN the RepEnrich2_subset.py to output discrete files for uniquely and multi-mapping reads
#python2.7 /work/GCRB/s157114/softwares/RepEnrich2/RepEnrich2_subset.py Sample.bam 30 Sample --pairedend TRUE


##Note LB Nature data is single-end. So modify above commands.
##First trying for WT H3K9me3 ChIP & Input only!!
gunzip -c /project/GCRB/Banaszynski_lab/shared/Yu-Ching/Repeat_Analysis/ChIP_data/LB_Nature/fastq_files/ESCH33WTH3K9me3.fastq.gz > ESCH33WTH3K9me3.fastq
bowtie2 -q -t -x /project/apps_database/iGenomes/Mus_musculus/UCSC/mm10/Sequence/Bowtie2Index/genome -U ESCH33WTH3K9me3.fastq -S ESCH33WTH3K9me3.sam
#Convert SAM to BAM
samtools view -bS ESCH33WTH3K9me3.sam > ESCH33WTH3K9me3.bam
#RUN the RepEnrich2_subset.py to output discrete files for uniquely and multi-mapping reads
python2.7 /work/GCRB/s157114/softwares/RepEnrich2/RepEnrich2_subset.py ESCH33WTH3K9me3.bam 30 Sample --pairedend FALSE

gunzip -c /project/GCRB/Banaszynski_lab/shared/Yu-Ching/Repeat_Analysis/ChIP_data/LB_Nature/fastq_files/ESCH33WTInput.fastq.gz > ESCH33WTInput.fastq
bowtie2 -q -t -x /project/apps_database/iGenomes/Mus_musculus/UCSC/mm10/Sequence/Bowtie2Index/genome -U ESCH33WTInput.fastq -S ESCH33WTInput.sam
#Convert SAM to BAM
samtools view -bS ESCH33WTInput.sam > ESCH33WTInput.bam
#RUN the RepEnrich2_subset.py to output discrete files for uniquely and multi-mapping reads
python2.7 /work/GCRB/s157114/softwares/RepEnrich2/RepEnrich2_subset.py ESCH33WTInput.bam 30 Sample --pairedend FALSE


##Step 4) Run RepEnrich2
###After we get the unique and multimap fastq, we run RepEnrich2 on the data
###Note 1: I have modified the original Repenrich2 script for compatibility with BioHPC
##Note 2: This step takes ages to run, so you might want to run one sample in each script (with 1 node and 240 hrs each).
python2.7 /work/GCRB/s157114/softwares/RepEnrich2/RepEnrich2_modified.py /work/GCRB/s157114/softwares/repeatmasker_seq/mm10_repeatmasker_without_simple_repeats.txt ./ESC_H33_WT_H3K9me3 ESC_H33_WT_H3K9me3 /work/GCRB/s157114/softwares/repeatmasker_seq/setup_folder_mm10_no_simple_repeats /project/GCRB/Banaszynski_lab/shared/Yu-Ching/Repeat_Analysis/ChIP_data/LB_Nature/RepEnrich_analysis/ESC_H33_WT_H3K9me3_R1_001_multimap.fastq /project/GCRB/Banaszynski_lab/shared/Yu-Ching/Repeat_Analysis/ChIP_data/LB_Nature/RepEnrich_analysis/ESC_H33_WT_H3K9me3_R1_001_unique.bam --pairedend FALSE --allcountmethod TRUE

##The output of Step 4 is a list of counts.txt files such as:
#SampleName_unique_counts.txt
#SampleName_total_counts.txt
#SampleName_fraction_counts.txt
#SampleName_family_total_counts.txt
#SampleName_family_fraction_counts.txt
#SampleName_class_total_counts.txt
#SampleName_class_fraction_counts.txt


##Step 5) Normalization of counts (club ChIPs in one group and Inputs in a separate group)
#Using specific counts files (either family or class depending on what you want), this step uses DESeq2 for normalization of raw counts from the fastq files.

#Load library (if you don’t have it in your R, then install it first)
library(“DESeq2”)

##ChIPs
A704_tSETD2 <- read.delim(‘A704_tSETD2_K36me3_family_fraction_counts.txt’, header=FALSE)
A704_Vec <- read.delim(‘A704_Vec_K36me3_family_fraction_counts.txt’, header=FALSE)

##Inputs
A704_tSETD2_Input<-read.delim(‘A704_tSETD2_K36me3_Input_family_fraction_counts.txt’, header=FALSE)
A704_Vec_Input <- read.delim(‘A704_Vec_K36me3_Input_family_fraction_counts.txt’, header=FALSE)
counts <- data.frame(row.names = A704_tSETD2[,1], A704_tSETD2 = round(A704_tSETD2[,2]), A704_Vec = round(A704_Vec[,2]))

colData<-data.frame(condition=c(“A704_tSETD2",“A704_Vec”))
row.names(colData)<-c(“A704_tSETD2",“A704_Vec”)
dds<-DESeqDataSetFromMatrix(countData=counts,colData=colData,design=~condition)
dds <- estimateSizeFactors(dds)
sizeFactors(dds)
dds<-dds[ rowSums(counts(dds)) > 1, ]
normalized_counts <- counts(dds, normalized=TRUE)
normalized_samples_fractions<-sweep(normalized_counts,MARGIN=2,FUN=“/”,STATS=colSums(normalized_counts))
normalized_counts_fractions<-do.call(“cbind”, list(normalized_counts,normalized_samples_fractions))


##Step 6) Getting fraction (as %) of a repeat element from the sum/pool of counts at all elements and then calculating enrichment over input.


##Step 7) Making representative figures such as heat maps.


##SETUP WITHOUT SIMPLE REPEATS
grep -v "Simple_repeat" mm10_repeatmasker.txt > mm10_repeatmasker_without_simple_repeats.txt
python2.7 /work/GCRB/s157114/softwares/RepEnrich2/RepEnrich2/RepEnrich2_setup.py ./mm10_repeatmasker_without_simple_repeats.txt /project/apps_database/iGenomes/Mus_musculus/UCSC/mm10/Sequence/Bowtie2Index/genome.fa ./setup_folder_mm10_no_simple_repeats/
