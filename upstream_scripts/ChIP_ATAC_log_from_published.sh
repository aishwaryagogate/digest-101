#!/bin/tcsh

#SBATCH --job-name=ChIP_ATAC_log                            # job name
#SBATCH --partition=super                                # select partion from 128GB, 256GB, 384GB, GPU and super
#SBATCH --nodes=1                                        # number of nodes requested by user
#SBATCH --time=0-72:00:00                                 # run time, format: D-H:M:S (max wallclock time)
#SBATCH --output=serialJob.%j.out                        # standard output file name
#SBATCH --error=serialJob.%j.time                        # standard error output file name
#SBATCH --mail-user=aishwarya.gogate@utsouthwestern.edu      # specify an email address
#SBATCH --mail-type=ALL                                  # send email when job status change (start, end, abortion and etc.)

set INDEX=/project/apps_database/iGenomes/Mus_musculus/UCSC/mm10/Sequence/BWAIndex/genome.fa
#set INDEX=/project/apps_database/iGenomes/Homo_sapiens/UCSC/hg19/Sequence/BWAIndex/genome.fa
set NUM_THREADS=30
set MISMATCH_PENALTY=8
set FILE=/project/GCRB/Banaszynski_lab/shared/Ryan1/ZFP_data/fastq_files

module load BWA/0.7.5
module load samtools
module load picard/1.117
module load bedtools
module load UCSC_userApps
module load macs/1.4.2
module load igvtools/2.3.71

mkdir bam_files

foreach sample(\
  file1\
  file2\
	)
    bwa mem\
       -t $NUM_THREADS\
       -B $MISMATCH_PENALTY\
       -M\
       $INDEX\
       $FILE/$sample.fastq.gz\
       | samtools view -q10 -bS -o - -\
       | samtools sort - -o ./bam_files/$sample.bam

    #samtools merge ./bam_files/$sample.bam\
	#$sample.L001.bam\
	#$sample.L002.bam\
	#$sample.L003.bam\
	#$sample.L004.bam

    samtools sort -m 8000000000 ./bam_files/$sample.bam

    #rm $sample.L001.bam $sample.L002.bam $sample.L003.bam $sample.L004.bam

    set MARK_DUP = /cm/shared/apps/picard/1.117/MarkDuplicates.jar

    mkdir nodup_files

    java -Xmx128g -jar $MARK_DUP\
    INPUT=./bam_files/$sample.bam\
    OUTPUT=./nodup_files/$sample.nodup.bam\
    METRICS_FILE=./nodup_files/metrics.$sample.txt\
    REMOVE_DUPLICATES=true\
    ASSUME_SORTED=true\
    TMP_DIR=temp_dir.$sample

    samtools index ./nodup_files/$sample.nodup.bam

    set CHROM_SIZE = /project/apps_database/iGenomes/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa.fai
#    set CHROM_SIZE = /project/apps_database/iGenomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa.fai

    mkdir bw_files

    bedtools genomecov -ibam ./nodup_files/$sample.nodup.bam\
			-bga\
			> ./bw_files/$sample.nodup.bedGraph

    sort -k1,1 -k2,2n ./bw_files/$sample.nodup.bedGraph > ./bw_files/$sample.nodup.sorted.bedGraph

    bedGraphToBigWig ./bw_files/$sample.nodup.sorted.bedGraph\
		     $CHROM_SIZE\
		     ./bw_files/$sample.nodup.bw
end

mkdir peaks_calls

    	macs14 -t ./nodup_files/file1.nodup.bam -g 1.87e9 -n ./peaks_calls/file1_peaks
      macs14 -t ./nodup_files/file2.nodup.bam -g 1.87e9 -n ./peaks_calls/file2_peaks
