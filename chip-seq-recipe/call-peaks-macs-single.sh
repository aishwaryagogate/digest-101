#!/bin/bash
# call-peaks-macs.sh

script_name="call-peaks-macs2.sh"
script_ver="2.1.0"

#Help function
usage() {
  echo "-h  --Help documentation for $script_name"
  echo "-f  --File path to experiment tagAlign files."
  echo "-s  --File path to control tagAlign files."
  echo "-x  --File path to experiment cross-correlation scores."
  echo "-r  --UCSC Reference genome (e.g. hg19, mm10)"
  echo "-o  --Path to output directory"
  echo "-v  --Version of script"
  echo "Example: $script_name -f 'foo1.tagAlign.gz' -s 'con1.tagAlign.gz' -x 'foo1.cc' -r 'hg19' [-o '/path/to/output/dir/']"
  exit 1
}

# Version function
version(){
    echo "$script_name $script_ver"
    exit 1
}


# Peak calling function
call_peak() {
    # Establish variables

    experiment=$1
    control=$2
    xcor_scores_input=$3
    chrom_sizes=$4
    genomesize=$5
    out_dir=$6

    #Extract the fragment length estimate from cross-correlation scores file
    fraglen=`cat $xcor_scores_input |  grep "predicted" | cut -f6 -d ' '`
    echo $fraglen

    # Generate narrow peaks and preliminary signal tracks
    base_fn=$(basename "${experiment}")
    prefix=${base_fn%.tagAlign.gz}
    macs2 callpeak -t $experiment -c $control -f BED -n $out_dir/$prefix -g $genomesize -p 1e-5 --nomodel --shift 0 --extsize $fraglen --keep-dup all -B --SPMR

    # Generate fold enrichment signal tracks
    macs2 bdgcmp -t $out_dir/$prefix\_treat_pileup.bdg -c $out_dir/$prefix\_control_lambda.bdg --outdir $out_dir -o $prefix\_FE.bdg -m FE

    # Genearte bigWigs from bedgraph to support vizualization
    bedtools slop -i $out_dir/$prefix\_FE.bdg -g /$chrom_sizes -b 0 | bedClip stdin $chrom_sizes $out_dir/$prefix.fc.signal.bedgraph
    bedSort $out_dir/$prefix.fc.signal.bedgraph $out_dir/$prefix.sorted.fc.signal.bedgraph
    bedGraphToBigWig $out_dir/$prefix.sorted.fc.signal.bedgraph $chrom_sizes $out_dir/$prefix.fc_signal.bw
    rm $out_dir/$prefix.fc.signal.bedgraph
    rm $out_dir/$prefix.sorted.fc.signal.bedgraph

    # Generate fold enrichment signal tracks

    # Compute sval = min(no. of reads in ChIP, no. of reads in control) / 1,000,000
    experiment_reads=`gzip -dc $experiment | wc -l | cut -f1`
    control_reads=`gzip -dc $control | wc -l | cut -f1`
    if [ $experiment_reads -ge $control_reads ]; then
        min=$control_reads
    else
        min=$experiment_reads
    fi
    sval=`echo $min/1000000 | bc -l`

    macs2 bdgcmp -t $out_dir/$prefix\_treat_pileup.bdg -c $out_dir/$prefix\_control_lambda.bdg --outdir $out_dir -o $prefix\_ppois.bdg -m ppois -S $sval

    # Genearte bigWigs from bedgraph to support vizualization
    bedtools slop -i $out_dir/$prefix\_ppois.bdg -g $chrom_sizes -b 0 | bedClip stdin $chrom_sizes $out_dir/$prefix.pval.signal.bedgraph
    bedSort $out_dir/$prefix.pval.signal.bedgraph $out_dir/$prefix.sorted.pval.signal.bedgraph
    bedGraphToBigWig $out_dir/$prefix.sorted.pval.signal.bedgraph $chrom_sizes $out_dir/$prefix.pval_signal.bw
    rm $out_dir/$prefix.pval.signal.bedgraph
    rm $out_dir/$prefix.sorted.pval.signal.bedgraph
}

main(){

    # Load required modules
    module load python/2.7.x-anaconda
    module load R/3.1.0-intel
    module load macs/2.1.0-20151222
    module load gcc/4.8.1
    module load bedtools/2.17.0
    module load UCSC_userApps/v317

    # Parsing options
    OPTIND=1 # Reset OPTIND
    while getopts :f:s:x:r:o:hv opt
        do
            case $opt in
                f) aln=$OPTARG;;
                s) control=$OPTARG;;
                x) aln_cross=$OPTARG;;
                r) ucsc_reference=$OPTARG;;
                o) out=$OPTARG;;
                h) usage;;
                v) version;;
            esac
        done

    shift $(($OPTIND -1))

    # Check for mandatory options
    if [[ -z $aln ]] || [[ -z $control ]] || [[ -z $ucsc_reference ]] || [[ -z $aln_cross ]]; then
        usage
    fi

    # Check if length of arguments in aln1, aln2 and exp are the same
    array_aln=(${aln//[,| ]/ })
    array_control=(${control//[,| ]/ })
    array_aln_cross=(${aln_cross//[,| ]/ })
    # Define the genome size to use
    if [ $ucsc_reference = 'hg19' ]; then
        genome_size='hs'
    elif [ $ucsc_reference = 'mm10' ]; then
        genome_size='mm'
    elif [ $ucsc_reference = 'mm9' ]; then
        genome_size='mm'
    else
        usage
    fi

    # Define the output directory, if none defined make the location relative to first file
    if [ -z $out ]; then
        one_parent=$(dirname "${aln}")
        out_dir=$(dirname "${one_parent}")\/$script_name-$script_ver/
    else
        out_dir=$out\/$script_name-$script_ver
    fi

    if [ ! -d $out_dir ]; then
        mkdir $out_dir
    fi

    # Align if file doesn't exist
    if [ ! -e $out_dir/metadata.json ]; then

        # split files
        rep1=${aln}
        con1=${control}
        rep1_xcor=${aln_cross}

        echo "* Downloading chrom.sizes..."
        fetchChromSizes $ucsc_reference > $out_dir/chrom.sizes
        chrom_sizes=$out_dir/chrom.sizes

        # Call peaks Replicates 1 and 2
        call_peak $rep1 $con1 $rep1_xcor $chrom_sizes $genome_size $out_dir

        # Remove chrom sizes
        rm $out_dir/chrom.sizes

        # Get input and output files and then print out metadata.json file
        input_files=("${array_aln[@]}" "${array_control[@]}" "${array_aln_cross[@]}" "$ucsc_reference")
        printf -v input "\"%s\"," "${input_files[@]}"
        input=${input%,}
        output_file=($out_dir\/*)
        printf -v output "\"%s\"," "${output_file[@]}"
        output=${output%,}
        printf '{"script name":"%s","script version":"%s", "input files": [%s], "output files": [%s]}' "$script_name" "$script_ver" "$input" "$output"  | python -m json.tool > $out_dir/metadata.json
    else
        echo "* Peaks have been called"
    fi
}

main "$@"

