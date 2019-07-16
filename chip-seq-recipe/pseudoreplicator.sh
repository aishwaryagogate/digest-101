#!/bin/bash
# pseudoreplicator.sh

script_name="pseudoreplicator.sh"
script_ver="1.0.0"

#Help function
usage() {
  echo "-h  --Help documentation for $script_name"
  echo "-f  --File path to filtered tagAlign"
  echo "-o  --Path to output directory"
  echo "-v  --Version of script"
  echo "Example: $script_name -f 'foo1.tagAlign' -o '/path/to/output/dir/'"
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
    module load R/3.2.1-intel
    module load UCSC_userApps/v317

    # Parsing options
    OPTIND=1 # Reset OPTIND
    while getopts :f:o:vh opt
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
    if [[ -z $aln_file ]] || [[ -z $out ]]; then
        usage
    fi


    # Define the output directory, if none defined make the location relative to first file
    if [ -z $out ]; then
        out_dir=$(dirname "${aln_file}")\/$script_name-$script_ver/
    else
        out_dir=$out\/$script_name-$script_ver
    fi

    if [ ! -d $out_dir ]; then
        mkdir $out_dir
    fi

    # Define the output file name, based on the first file
    aln_fn=$(basename "${aln_file}")
    output_fp=${aln_fn%.tagAlign.gz}

    if [ ! -d $out_dir/$output_fp\.metadata.json ]; then

        # Count lines
        total_lines=`gzip -dc $aln_file | wc -l`

        #number of lines in each split
        sum=$(expr $total_lines + 1)
        nlines=$(expr $sum / 2)

        # Split files
        splits_prefix='temp_split'
        gzip -dc $aln_file | shuf | split -a 2 -d -l $nlines - $out_dir/$splits_prefix
        gzip -c $out_dir/$splits_prefix\00 > $out_dir/$output_fp.pr1.tagAlign.gz
        gzip -c $out_dir/$splits_prefix\01 > $out_dir/$output_fp.pr2.tagAlign.gz
        rm $out_dir/$splits_prefix\00
        rm $out_dir/$splits_prefix\01

        # Get input and output files and then print out metadata.json file
        input_files=("$aln_file")
        printf -v input "\"%s\"," "${input_files[@]}"
        input=${input%,}
        output_file=($out_dir\/$output_fp*)
        printf -v output "\"%s\"," "${output_file[@]}"
        output=${output%,}
        printf '{"script name":"%s","script version":"%s", "input files": [%s], "output files": [%s]}' "$script_name" "$script_ver" "$input" "$output"  | python -m json.tool > $out_dir/$output_fp.metadata.json
        echo "* Finished."
    else
        echo "* $aln_fn pseudoreplication has been run. "
    fi
}

main "$@"
