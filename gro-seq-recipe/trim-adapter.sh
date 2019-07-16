#!/bin/bash
# trim-adapter.sh

script_name="trim-adapter.sh"
script_ver="1.0.0"

#Help function
usage() {
  echo "-h  --Help documentation for $script_name"
  echo "-f  --Path to fastq file"
  echo "-a  --3' adapter sequence (e.g. TCGTATCCCGTCTTCTGCTTG)"
  echo "-p  --PolyA sequence (default AAAAAAAAAAAAAAAAAAAA)"
  echo "-o  --Path to output directory"
  echo "-v  --Version of script"
  echo "Example: $script_name -f 'foo1.fastq.gz' -a 'TCGTATCCCGTCTTCTGCTTG' [-p 'AAAAAAAAAAA'] [-o '/path/to/output/dir/']"
  exit 1
}

# Version function
version(){
    echo "$script_name $script_ver"
    exit 1
}

main(){

    # Load required modules
    module load cutadapt/1.9.1

    # Parsing options
    OPTIND=1 # Reset OPTIND
    while getopts :f:a:povh opt
        do
            case $opt in
                f) raw_file=$OPTARG;;
                a) adapter=$OPTARG;;
                p) polya=$OPTARG;;
                o) out=$OPTARG;;
                v) version;;
                h) usage;;
            esac
        done

    shift $(($OPTIND -1))

    # Check for mandatory options
    if [[ -z $raw_file ]] || [[ -z $adapter ]]; then
        usage
    fi

    # Define PolyA sequence to trim
    if [[ -z $polya ]]; then
        polya='AAAAAAAAAAAAAAAAAAAA'
    fi

    # Define the output directory
    raw_fn=$(basename "${raw_file}")
    raw_prefix=${raw_fn%.fastq.gz}
    if [ -z $out ]; then
        out_dir=$(dirname "${raw_file}")\/$script_name-$script_ver
    else
        out_dir=$out\/$script_name-$script_ver
    fi

    # Make Out directory if doesn't exist
    if [ ! -e $out_dir ]; then
        mkdir $out_dir
    fi

    # Run cutadapt
    if [ ! -e $out_dir\/$raw_prefix.noIllAdapt.fastq.gz ]; then
        echo "* Triming 3' adapter sequence '$raw_fn'"
        cutadapt -a $adapter -z -e 0.10 --minimum-length=32 --output=$out_dir/$raw_prefix.noIllAdapt.fastq.gz $raw_file

        echo "* Triming PolyA sequence '$raw_fn'"
        cutadapt -a $polya -z -e 0.10 --minimum-length=32 --output=$out_dir/$raw_prefix.noAdapt.fastq.gz $out_dir/$raw_prefix.noIllAdapt.fastq.gz


        # Get input and output files and then print out metadata.json file
        input_files=("$raw_file")
        printf -v input "\"%s\"," "${input_files[@]}"
        input=${input%,}
        output_file=($out_dir\/$raw_prefix*)
        printf -v output "\"%s\"," "${output_file[@]}"
        output=${output%,}
        printf '{"script name":"%s","script version":"%s", "input files": [%s], "output files": [%s]}' "$script_name" "$script_ver" "$input" "$output"  | python -m json.tool > $out_dir/$raw_fn.metadata.json
        echo "* Finished."
    else
        echo "* Adapters have been trimmed '$raw_fn'"
    fi
}

main "$@"
