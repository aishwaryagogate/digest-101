#!/bin/tcsh

#SBATCH --job-name=fingerprint                            # job name
#SBATCH --partition=super                                # select partion from 128GB, 256GB, 384GB, GPU and super
#SBATCH --nodes=1                                        # number of nodes requested by user
#SBATCH --time=0-24:00:00                                 # run time, format: D-H:M:S (max wallclock time)
#SBATCH --output=serialJob.%j.out                        # standard output file name
#SBATCH --error=serialJob.%j.time                        # standard error output file name
#SBATCH --mail-user=aishwarya.gogate@utsouthwestern.edu      # specify an email address
#SBATCH --mail-type=ALL                                  # send email when job status change (start, end, abortion and etc.)

module load samtools
module load bedtools
module load deeptools

##Some useful info:
##Can be used for QC of data.
##This command determines how well the signal in the ChIP-seq sample can be differentiated from the background distribution of reads in the control sample.
##For factors that will enrich well-defined, rather narrow regions (e.g. transcription factors such as p300), the resulting plot can be used to assess the strength of a ChIP.
##But the broader the enrichments are to be expected, the less clear the plot will be. ALso, if you do not know what kind of signal to expect, this plot is helpful.
##It will give you a straight-forward indication of how careful you will have to be during your downstream analyses to separate biological noise from meaningful signal.

##What the plot tells you?
##An ideal [input] with perfect uniform distribution of reads along the genome (i.e. without enrichments in open chromatin etc.) should generate a straight diagonal line.
##A very specific and strong ChIP enrichment will be indicated by a prominent and steep rise of the cumulative sum towards the highest rank.
##This means that a big chunk of reads from the ChIP sample is located in few bins which corresponds to high, narrow enrichments typically seen for TFs.

#Example usage: plotFingerprint -b treatment.bam control.bam -plot fingerprint.pdf
#See here for optional arguments: https://deeptools.readthedocs.io/en/2.1.0/content/tools/plotFingerprint.html#what-the-plots-tell-you
#Note: BAM files need to be indexed!
plotFingerprint -b 129_H3K27ac_SI_S1_R1.sorted.bam 129_INPUT_SI_S5_R1.sorted.bam minMappingQuality 30 --skipZeros --plotTitle= "Fingerprints of WT & Input H3K27ac" --plotFile WT_fingerprints.pdf --outRawCounts WT_fingerprints.tab
