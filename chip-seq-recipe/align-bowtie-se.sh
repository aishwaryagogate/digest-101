#!/bin/bash
# align-bowtie-se.sh

script_name="align-bowtie-se.sh"
script_ver="1.0.0"

#Help function
usage() {
  echo "-h Help documentation for $script_name"
  echo "-f  --Path to fastq file"
  echo "-r  --UCSC Reference genome (e.g. hg19, mm10)"
  echo "-o  --Path to output directory"
  echo "-v  --Version of script"
  echo "Example: $script_name -f 'foo_1.fastq.gz' -r 'hg19' [-o '/path/to/output/dir/']"
  exit 1
}

# Version function
version(){
    echo "$script_name $script_ver"
    exit 1
}

main(){

    # Load required modules
    module load bowtie/1.0.0
    module load samtools/0.1.19
    module load bedtools/2.17.0
    module load iGenomes/2013-03-25

    # Parsing options
    OPTIND=1 # Reset OPTIND
    while getopts :f:r:o:vh opt
        do
            case $opt in
                f) raw_file=$OPTARG;;
                r) ucsc_reference=$OPTARG;;
                o) out=$OPTARG;;
                v) version;;
                h) usage;;
            esac
        done

    shift $(($OPTIND -1))

    # Check for mandatory options
    if [[ -z $raw_file ]] || [[ -z $ucsc_reference ]]; then
        usage
    fi

    # Define the index to use
    # Establish the genome file version to use and BOWTIEINDEX
    if [ $ucsc_reference = 'hg19' ]; then
        index=$iGENOMES_DB_DIR\/Homo_sapiens/UCSC/$ucsc_reference/Sequence/BowtieIndex/genome
    elif [ $ucsc_reference = 'mm10' ]; then
        index=$iGENOMES_DB_DIR\/Mus_musculus/UCSC/$ucsc_reference/Sequence/BowtieIndex/genome
    elif [ $ucsc_reference = 'mm9' ]; then
        index=$iGENOMES_DB_DIR\/Mus_musculus/UCSC/$ucsc_reference/Sequence/BowtieIndex/genome
    else
        usage
    fi

    # Define the output directory, if none defined make the location relative to first file
    if [ -z $out ]; then
        out_dir=$(dirname "${raw_file}")\/$script_name-$script_ver/
    else
        out_dir=$out\/$script_name-$script_ver
    fi

    if [ ! -d $out_dir ]; then
        mkdir $out_dir
    fi

    # Define the output file name, based on the first file
    raw_fn=$(basename "${raw_file}")
    output_fp=${raw_fn%_*}

    # Align if file doesn't exist
    if [ ! -e $out_dir/$output_fp.sorted.bam ]; then

        #unizip files
        file_un_zip=${raw_fn%.gz}
        gunzip -c $raw_file > $out_dir/$file_un_zip

        bowtie -n 2 -k 1 -m 1 -S -t $index $out_dir/$file_un_zip $out_dir/$output_fp.sam

        # Convert sam to bam
        samtools view -bh -S $out_dir/$output_fp.sam > $out_dir/$output_fp.unsorted.bam

        # Sort bam
    	samtools sort $out_dir/$output_fp.unsorted.bam $out_dir/$output_fp.sorted

    	rm $out_dir/$output_fp.unsorted.bam
        rm $out_dir/$output_fp.sam
        rm $out_dir/$file_un_zip

        # Get input and output files and then print out metadata.json file
        input_files=("$index" "$raw_file")
        printf -v input "\"%s\"," "${input_files[@]}"
        input=${input%,}
        output_file=($out_dir\/$output_fp*)
        printf -v output "\"%s\"," "${output_file[@]}"
        output=${output%,}
        printf '{"script name":"%s","script version":"%s", "input files": [%s], "output files": [%s]}' "$script_name" "$script_ver" "$input" "$output"  | python -m json.tool > $out_dir/$output_fp.metadata.json

    else
        echo "* $output_fp has already been aligned"
    fi
}

main "$@"
