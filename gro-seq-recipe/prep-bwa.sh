#!/bin/bash
# prep-bwa.sh

script_name="prep-bwa.sh"
script_ver="1.0.0"

#Help function
usage() {
  echo "-h  --Help documentation for $script_name"
  echo "-r  --Path to reference genome"
  echo "-o  --Output file name"
  echo "-v  --Version of script"
  echo "Example: $script_name -r '/path/to/reference/genome.fa' -o 'output_file'"
  exit 1
}

# Version function
version(){
    echo "$script_name $script_ver"
    exit 1
}


main(){

    # Load required modules
    module load bwa/intel/0.7.12
    module load samtools/0.1.19
    module load picard/1.117

    # Parsing options
    OPTIND=1 # Reset OPTIND
    while getopts :r:o:vh opt
        do
            case $opt in
                r) reference_file=$OPTARG;;
                o) output_file=$OPTARG;;
                v) version;;
                h) usage;;
            esac
        done

    shift $(($OPTIND -1))

    # Check for mandatory options
    if [[ -z $reference_file ]] || [[ -z $output_file ]]; then
        usage
    fi

    reference_fn=$(basename "${reference_file}")
    parent=$(dirname "${reference_file}")
    out_dir=$(dirname "${parent}")\/$script_name-$script_ver/
    if [ ! -e $out_dir/$output_file ]; then
        mkdir $out_dir
        cp $reference_file $out_dir/

        echo "* Value of reference genome: '$reference_file'"
        # Generate the BWA index
        echo "* Generating bwa index genome."
        bwa index -a bwtsw $out_dir/$reference_fn

        # Generate the fasta file index
        echo "* Generate fasta file index."
        samtools faidx $out_dir/$reference_fn

        # Generate the sequence index
        echo "* Generate sequence index."
        java -jar CreateSequenceDictionary.jar REFERENCE=$out_dir/$reference_fn OUTPUT=$out_dir/$output_file

        # Get input and output files and then print out metadata.json file
        input_files=("$reference_file")
        printf -v input "\"%s\"," "${input_files[@]}"
        input=${input%,}
        output_files=($out_dir/*)
        printf -v output "\"%s\"," "${output_file[@]}"
        output=${output%,}
        printf '{"script name":"%s","script version":"%s", "input files": [%s], "output files": [%s]}' "$script_name" "$script_ver" "$input" "$output"  | python -m json.tool > $out_dir/$alignment_fn.metadata.json
        echo "* Finished."
    else
        echo "* Index has been generated"
    fi
}

main "$@"
