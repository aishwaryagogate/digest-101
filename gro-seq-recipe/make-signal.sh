#!/bin/bash
# make-signal.sh

script_name="make-signal.sh"
script_ver="1.0.0"

#Help function
usage() {
  echo "-h Help documentation for $script_name"
  echo "-f  --File path to Replicate 1 alignments"
  echo "-s  --File path to Replicate 2 alignments"
  echo "-r  --UCSC Reference genome (e.g. hg19, mm10)"
  echo "-p  --The prefix of the signal files"
  echo "-o  --Path to output directory"
  echo "-v  --Version of script"
  echo "Example: $script_name -f 'foo_1.bam' [-s 'foo_2.bam'] -r 'hg19' -p 'foo' -o '/path/to/output/dir/'"
  exit 1
}

# Version function
version(){
    echo "$script_name $script_ver"
    exit 1
}


main(){

    # Load required modules
    module load python/2.7.x-anaconda
    module load R/3.1.0-intel

    # Parsing options
    OPTIND=1 # Reset OPTIND
    while getopts :f:s:r:p:o:vh opt
        do
            case $opt in
                f) aln1=$OPTARG;;
                s) aln2=$OPTARG;;
                r) ucsc_reference=$OPTARG;;
                p) prefix=$OPTARG;;
                o) out=$OPTARG;;
                v) version;;
                h) usage;;
            esac
        done

    shift $(($OPTIND -1))

    # Check for mandatory options
    if [[ -z $aln1 ]] || [[ -z $ucsc_reference ]] || [[ -z $prefix ]] || [[ -z $out ]]; then
        usage
    fi

    # Define the out directory
    out_dir=$out\/$script_name-$script_ver

    # Make sure directories exist
    if [ ! -e $out ]; then
        mkdir $out
    fi

    if [ ! -e $out_dir ]; then
        mkdir $out_dir
    fi

    # Run make-signal.R if files doesn't exist
    if [ ! -e $out_dir/$prefix\_Plus.bigWig ]; then
        # Check to see if 2 replicates
        if [[ -z $aln2 ]]; then
            Rscript make-signal.R --replicate1 $aln1 --genome $ucsc_reference --prefix $prefix --out $out_dir

            # Set input file
            input_files=("$aln1")
        else
            Rscript make-signal.R --replicate1 $aln1 --replicate2 $aln2 --genome $ucsc_reference --prefix $prefix --out $out_dir

            # Set input files
            input_files=("$aln1" "$aln2")
        fi

        # Get input and output files and then print out metadata.json file
        printf -v input "\"%s\"," "${input_files[@]}"
        input=${input%,}
        output_file=($out_dir\/$prefix*)
        printf -v output "\"%s\"," "${output_file[@]}"
        output=${output%,}
        printf '{"script name":"%s","script version":"%s", "input files": [%s], "output files": [%s]}' "$script_name" "$script_ver" "$input" "$output"  | python -m json.tool > $out_dir/$prefix\_metadata.json
    else
        if [[ -z $aln2 ]]; then
            aln1_fn=$(basename "$aln1")
            echo "* Signal has been made from $aln1_fn "
        else
            aln1_fn=$(basename "$aln1")
            aln2_fn=$(basename "$aln2")
            echo "* Signal has been made from $aln1_fn and $aln2_fn "
        fi
    fi
}

main "$@"
