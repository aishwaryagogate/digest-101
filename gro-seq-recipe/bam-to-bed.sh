#!/bin/bash
# bam-to-bed.sh

script_name="bam-to-bed.sh"
script_ver="1.0.0"

#Help function
usage() {
  echo "-h  --Help documentation for $script_name"
  echo "-f  --Path to bam file"
  echo "-o  --Path to output directory"
  echo "-v  --Version of script"
  echo "Example: $script_name -f 'foo1.bam' [-o '/path/to/output/dir/']"
  exit 1
}

# Version function
version(){
    echo "$script_name $script_ver"
    exit 1
}

main(){

    # Load required modules
    module load bedtools/2.17.0

    # Parsing options
    OPTIND=1 # Reset OPTIND
    while getopts :f:ovh opt
        do
            case $opt in
                f) aln_file=$OPTARG;;
                o) out=$OPTARG;;
                v) version;;
                h) usage;;
            esac
        done

    shift $(($OPTIND -1))

    # Check for mandatory options
    if [[ -z $aln_file ]]; then
        usage
    fi

    # Define the output directory
    aln_fn=$(basename "${aln_file}")
    aln_prefix=${raw_fn%.bam}
    if [ -z $out ]; then
        out_dir=$(dirname "${aln_file}")\/$script_name-$script_ver
    else
        out_dir=$out\/$script_name-$script_ver
    fi

    # Make Out directory if doesn't exist
    if [ ! -e $out_dir ]; then
        mkdir $out_dir
    fi

        # Run bam to bed conversion
        if [ ! -e $out_dir\/$aln_prefix.bed ]; then
            echo "* Converting '$aln_fn' to bed"
            bamTobed -i $aln_file > $out_dir/$aln_prefix.bed

            # Get input and output files and then print out metadata.json file
            input_files=("$aln_file")
            printf -v input "\"%s\"," "${input_files[@]}"
            input=${input%,}
            output_file=($out_dir\/$aln_prefix*)
            printf -v output "\"%s\"," "${output_file[@]}"
            output=${output%,}
            printf '{"script name":"%s","script version":"%s", "input files": [%s], "output files": [%s]}' "$script_name" "$script_ver" "$input" "$output"  | python -m json.tool > $out_dir/$aln_prefix.metadata.json
            echo "* Finished."
        else
            echo "* BED file has been generated for '$aln_fn'"
        fi
}

main "$@"
