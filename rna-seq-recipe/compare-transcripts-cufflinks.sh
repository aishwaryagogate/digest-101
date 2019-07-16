#!/bin/bash
# compare-transcripts-cufflinks.sh

script_name="compare-transcripts-cufflinks.sh"
script_ver="1.0.0"

#Help function
usage() {
  echo "-h  --Help documentation for $script_name"
  echo "-g  --File path to transcript files, space seperated"
  echo "-a  --File path to reference annotation"
  echo "-o  --Path to the output directory"
  echo "-v  --Version of script"
  echo "Example: $script_name -g '/path/to/transcripts/foo_rep1.gtf /path/to/transcripts/foo_rep2.gtf' -a '/path/to/gtf/gencode.v19.annotation.gtf' -o 'experiment_152_comparison'"
  exit 1
}

# Version function
version(){
    echo "$script_name $script_ver"
    exit 1
}

main(){

    # Load required modules
    module load cufflinks/2.2.1
    module load python/2.7.x-anaconda

    # Parsing options
    OPTIND=1 # Reset OPTIND
    while getopts :g:a:o:vh opt
        do
            case $opt in
                g) transcripts=$OPTARG;;
                a) annotation=$OPTARG;;
                o) out=$OPTARG;;
                v) version;;
                h) usage;;
            esac
        done

    shift $(($OPTIND -1))

    # Check for mandatory options
    if [[ -z $transcripts ]] || [[ -z $annotation ]] || [[ -z $out ]]; then
        usage
    fi

    # Make Out directory if doesn't exist
    if [ ! -e $out ]; then
        mkdir $out
    fi

    # Compare transcripts if not already done so
    current_directory=$(pwd)
    out_dir=$out\/$script_name-$script_ver
    if [ ! -d $out_dir ]; then
        mkdir $out_dir

        echo "* Value of annotation: '$annotation'"
        echo "* Generating Comparison... "
        cuffcompare -p 10 -r $annotation ${transcripts}

        # Rename and move files
        echo "* Moving Files... "
        rename_files=("cuffcmp.combined.gtf" "cuffcmp.tracking" "cuffcmp.stats" "cuffcmp.loci")
        for file in ${rename_files[@]}
        do
            mv $current_directory\/$file $out_dir\/$file
        done

        for file in ${transcripts[@]}
        do
            file_fn=$(basename "${file}")
            file_fp=$(dirname "${file}")
            mv  $file_fp/cuffcmp.$file_fn.refmap  $file_fp/cuffcmp.$file_fn.tmap $out_dir\/
        done


        # Get input and output files and then print out metadata.json file
        array_transcripts=($transcripts)
        input_files=("$annotation" "${array_transcripts[@]}")
        printf -v input "\"%s\"," "${input_files[@]}"
        input=${input%,}
        output_file=($out_dir\/*)
        printf -v output "\"%s\"," "${output_file[@]}"
        output=${output%,}
        printf '{"script name":"%s","script version":"%s", "input files": [%s], "output files": [%s]}' "$script_name" "$script_ver" "$input" "$output"  | python -m json.tool > $out_dir/metadata.json
        echo "* Finished."
    else
        echo "* Transcripts comparisons have been made."
    fi
}

main "$@"
