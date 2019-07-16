#!/bin/bash
# call-combined-peaks-macs2.sh

script_name="call-combined-peaks-macs2"
script_ver="1.0.0"

#Help function
usage() {
  echo "-h Help documentation for $script_name"\\n
  echo "-a --Path to the alignment."
  echo "-c --Path to the controls."
  echo "-p --Prefix to use"
  echo "Example: $script_name -a 'foo_rep1.bed foo_rep2.bed' -c 'control_rep1.bed control_rep2.bed' -p 'foo'"
  exit 1
}

main(){

    #TODO: Refactor code

    # Load required modules
    module load python/2.7.x-anaconda
    module load gcc/4.8.1
    module load bedtools/2.17.0
    module load UCSC_userApps/v317

    # Source macs2
    source activate macs2

    # Parsing options
    OPTIND=1 # Reset OPTIND
    while getopts :a:c:p:h opt
        do
            case $opt in
                a)  BEDFILES=$OPTARG;;
                c)  INPUTBEDFILES=$OPTARG;;
                p)  PRE=$OPTARG;;
                h)  usage;;
            esac
        done

    shift $(($OPTIND -1))

    array_control=(${INPUTBEDFILES//[,| ]/ })
    array_files=(${BEDFILES//[,| ]/ })

    BEDFILE=${array_files[1]}
    BASEFILE=$(basename "${BEDFILE}")
    ONEPARENT=$(dirname "${BEDFILE}") # Get the first parent directory and go up
    OUTDIR=$(dirname "${ONEPARENT}")\/$script_name-$script_ver

    # Make Out directory if doesn't exist
    if [ ! -e $OUTDIR ]; then
        mkdir $OUTDIR
    fi

    # Align if file doesn't exist
    if [ ! -e $OUTDIR/$PRE\_peaks.xls ]; then

        # Generate narrow peaks
        macs2 callpeak -t $BEDFILES -c $INPUTBEDFILES -f BED -n $OUTDIR/$PRE -g hs -B --SPMR

        # Rescale Col5 scores to range 10-1000 to conform to narrowPeak.as format (score must be <1000)
        python rescale.py $OUTDIR/$PRE\_peaks.narrowPeak

        # Sort by Col8 in descending order and replace long peak names in Column 4 with Peak_<peakRank>
        sort -k 8gr,8gr $OUTDIR/rescaled-$PRE\_peaks.narrowPeak | awk 'BEGIN{OFS="\t"}{$4="Peak_"NR ; print $0}'| tee $OUTDIR/$PRE.narrowPeak
        sort -k1,1 -k2,2n  $OUTDIR/$PRE.narrowPeak >  $OUTDIR/$PRE.sorted.narrowPeak

        # Generate bigBeds from beds
        bedToBigBed $OUTDIR/$PRE.sorted.narrowPeak -type=bed6+4 -as=/home2/s163035/encValData/as/narrowPeak.as /home2/s163035/encValData/hg19/chrom.sizes  $OUTDIR/$PRE.narrowPeak.bb
        rm $OUTDIR/sorted-$PRE\_peaks.narrowPeak
        rm $OUTDIR/rescaled-$PRE\_peaks.narrowPeak
        rm $OUTDIR/$PRE.sorted.narrowPeak


        # Generate broad and gapped peaks
        macs2 callpeak -t $BEDFILE -c $INPUTBEDFILES -f BED -n $OUTDIR/$PRE -g hs --broad

        # Rescale Col5 scores to range 10-1000 to conform to broadPeak.as and gappedPeak.as format (score must be <1000)
        python rescale.py $OUTDIR/$PRE\_peaks.broadPeak
        python rescale.py $OUTDIR/$PRE\_peaks.gappedPeak

        # Sort by Col8 (for broadPeak) or Col 14(for gappedPeak)  in descending order and replace long peak names in Column 4 with Peak_<peakRank>
        sort -k 8gr,8gr $OUTDIR/rescaled-$PRE\_peaks.broadPeak | awk 'BEGIN{OFS="\t"}{$4="Peak_"NR ; print $0}'| tee $OUTDIR/$PRE.broadPeak
        sort -k1,1 -k2,2n  $OUTDIR/$PRE.broadPeak >  $OUTDIR/$PRE.sorted.broadPeak
        sort -k 14gr,14gr $OUTDIR/rescaled-$PRE\_peaks.gappedPeak | awk 'BEGIN{OFS="\t"}{$4="Peak_"NR ; print $0}'| tee $OUTDIR/$PRE.gappedPeak
        sort -k1,1 -k2,2n  $OUTDIR/$PRE.gappedPeak >  $OUTDIR/$PRE.sorted.gappedPeak

        # Generate bigBeds from beds
        bedToBigBed $OUTDIR/$PRE.sorted.broadPeak -type=bed6+3 -as=/home2/s163035/encValData/as/broadPeak.as /home2/s163035/encValData/hg19/chrom.sizes  $OUTDIR/$PRE.broadPeak.bb
        bedToBigBed $OUTDIR/$PRE.sorted.gappedPeak -type=bed12+3 -as=/home2/s163035/encValData/as/gappedPeak.as /home2/s163035/encValData/hg19/chrom.sizes  $OUTDIR/$PRE.gappedPeak.bb
        rm $OUTDIR/sorted-$PRE\_peaks.broadPeak
        rm $OUTDIR/rescaled-$PRE\_peaks.broadPeak
        rm $OUTDIR/$PRE.sorted.broadPeak
        rm $OUTDIR/sorted-$PRE\_peaks.gappedPeak
        rm $OUTDIR/rescaled-$PRE\_peaks.gappedPeak
        rm $OUTDIR/$PRE.sorted.gappedPeak


        # Generate fold enrichment signal tracks
        macs2 bdgcmp -t $OUTDIR/$PRE\_treat_pileup.bdg -c $OUTDIR/$PRE\_control_lambda.bdg --outdir $OUTDIR -o $PRE\_FE.bdg -m FE

        # Genearte bigWigs from bedgraph to support vizualization
        bedtools slop -i $OUTDIR/$PRE\_FE.bdg -g /home2/s163035/encValData/hg19/chrom.sizes -b 0 | bedClip stdin /home2/s163035/encValData/hg19/chrom.sizes $OUTDIR/$PRE.fc.signal.bedgraph
        bedSort $OUTDIR/$PRE.fc.signal.bedgraph $OUTDIR/$PRE.sorted.fc.signal.bedgraph
        bedGraphToBigWig $OUTDIR/$PRE.sorted.fc.signal.bedgraph /home2/s163035/encValData/hg19/chrom.sizes $OUTDIR/$PRE.fc_signal.bw
        rm $OUTDIR/$PRE.fc.signal.bedgraph
        rm $OUTDIR/$PRE.sorted.fc.signal.bedgraph

        # Get input and output files and then print out metadata.json file
        input_files=("${array_files[@]}" "${array_control[@]}")
        printf -v input "\"%s\"," "${input_files[@]}"
        input=${input%,}
        output_file=($OUTDIR\/$PRE*)
        printf -v output "\"%s\"," "${output_file[@]}"
        output=${output%,}
        printf '{"script name":"%s","script version":"%s", "input files": [%s], "output files": [%s]}' "$script_name" "$script_ver" "$input" "$output"  | python -m json.tool > $OUTDIR/$PRE.metadata.json

    else
        echo "* Mearged peaks have already been called"
    fi
}
