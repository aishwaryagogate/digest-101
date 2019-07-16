#!/bin/bash
# make-normal-signal.sh

script_name="make-normal-signal.sh"
script_ver="1.0.0"

#Help function
usage() {
  echo "-h Help documentation for $script_name"\\n
  echo "-f --File path to Replicate 1 alignments."
  echo "-s --File path to Replicate 2 alignments."
  echo "-r  --UCSC Reference genome (e.g. hg19, mm10)"
  echo "-e  --The experiment names."
  echo "-o  --The ouput directory."
  echo "Example: $script_name -f 'foo1.bam man1.bam' -s 'foo2.bam man2.bam' -r 'hg19' -e 'foo man' -o output"
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
    while getopts :f:s:r:e:o:hv opt
        do
            case $opt in
                f) aln1=$OPTARG;;
                s) aln2=$OPTARG;;
                r) ucsc_reference=$OPTARG;;
                e) exp=$OPTARG;;
                o) out=$OPTARG;;
                h) usage;;
                v) version;;
            esac
        done

    shift $(($OPTIND -1))

    # Check for mandatory options
    if [[ -z $aln1 ]] || [[ -z $aln2 ]] || [[ -z $ucsc_reference ]] || [[ -z $exp ]] || [[ -z $out ]]; then
        usage
    fi

    # Check if length of arguments in aln1, aln2 and exp are the same
    array_aln1=(${aln1//[,| ]/ })
    array_aln2=(${aln2//[,| ]/ })
    array_exp=(${aln2//[,| ]/ })
    if [[(${#array_aln1[@]} -ne ${#array_aln2[@]}) && (${#array_aln1[@]} -ne ${#array_exp[@]}) && (${#array_aln2[@]} -ne ${#array_exp[@]}) ]]; then
        echo "The number of arguments are not equal."
        exit 1
    fi

    # Define the out directory
    out_dir=$out\/$script_name-$script_ver

    # Make sure directories exist
    if [ ! -e $out ]; then
        mkdir $out
    fi

    # Run make-normal-signal.R if directory doesn't exist
    if [ ! -d $out_dir ]; then
        mkdir $out_dir
        Rscript make-normal-signal.R --replicate1 $aln1 --replicate2 $aln2 --experiment $exp --genome $ucsc_reference --out $out_dir

        # Get input and output files and then print out metadata.json file
        input_files=("${array_aln1[@]}" "${array_aln2[@]}")
        printf -v input "\"%s\"," "${input_files[@]}"
        input=${input%,}
        output_file=($out_dir\/*)
        printf -v output "\"%s\"," "${output_file[@]}"
        output=${output%,}
        printf '{"script name":"%s","script version":"%s", "input files": [%s], "output files": [%s]}' "$script_name" "$script_ver" "$input" "$output"  | python -m json.tool > $out_dir/metadata.json
    else
        echo "* Normalized Signal has been made. "
    fi
}

main "$@"
