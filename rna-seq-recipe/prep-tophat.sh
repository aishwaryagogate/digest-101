#!/bin/bash
# prep-tophat.sh

script_name="prep-tophat.sh"
script_ver="1.0.0"

#Help function
usage() {
  echo "-h  --Help documentation for $script_name"
  echo "-r  --UCSC Reference genome (e.g. hg19, mm10)"
  echo "-a  --Path to the annotation file"
  echo "-s  --Path to the Spike-in set (e.g. ERCC)"
  echo "-v  --Version of script"
  echo "Example: $script_name -r 'hg19' -a '/path/to/gtf/gencode.v19.annotation.gtf' [-s '/path/to/spike-in/ercc_92.fa']"
  exit 1
}

# Version function
version(){
    echo "$script_name $script_ver"
    exit 1
}

main(){

    # Load required modules
    module load perl/5.18.2
    module load bowtie2/intel/2.1.0
    module load tophat/2.0.12
    module load iGenomes/2013-03-25
    module load python/2.7.x-anaconda

    # Parsing options
    OPTIND=1 # Reset OPTIND
    while getopts :r:a:svh opt
        do
            case $opt in
                r) ucsc_reference=$OPTARG;;
                a) annotation=$OPTARG;;
                s) spike_in=$OPTARG;;
                v) version;;
                h) usage;;
            esac
        done

    shift $(($OPTIND -1))

    # Check for mandatory options
    if [[ -z $ucsc_reference ]] || [[ -z $annotation ]]; then
        usage
    fi

    # Establish the genome file version to use and BOWTIEINDEX
    if [ $ucsc_reference = 'hg19' ]; then
        bowtie2_index=$iGENOMES_DB_DIR\/Homo_sapiens/UCSC/$ucsc_reference/Sequence/Bowtie2Index/genome
        ucsc_genome_fn=$iGENOMES_DB_DIR\/Homo_sapiens/UCSC/$ucsc_reference/Sequence/WholeGenomeFasta/genome.fa
    elif [ $ucsc_reference = 'mm10' ]; then
        bowtie2_index=$iGENOMES_DB_DIR\/Mus_musculus/UCSC/$ucsc_reference/Sequence/Bowtie2Index/genome
        ucsc_genome_fn=$iGENOMES_DB_DIR\/Mus_musculus/UCSC/$ucsc_reference/Sequence/WholeGenomeFasta/genome.fa
    elif [ $ucsc_reference = 'mm9' ]; then
        bowtie2_index=$iGENOMES_DB_DIR\/Mus_musculus/UCSC/$ucsc_reference/Sequence/Bowtie2Index/genome
        ucsc_genome_fn=$iGENOMES_DB_DIR\/Mus_musculus/UCSC/$ucsc_reference/Sequence/WholeGenomeFasta/genome.fa
    elif [ $ucsc_reference = 'dm3' ]; then
        bowtie2_index=$iGENOMES_DB_DIR\/Drosophila_melanogaster/UCSC/$ucsc_reference/Sequence/Bowtie2Index/genome
        ucsc_genome_fn=$iGENOMES_DB_DIR\/Drosophila_melanogaster/UCSC/$ucsc_reference/Sequence/WholeGenomeFasta/genome.fa
    else
        usage
    fi

    # Establish the output directory relative to the
    annotation_fp=$(dirname "${annotation}")
    out_fn=$(basename "${annotation}")
    if grep -q merge-annotation <<< $annotation_fp; then
        out_dir=$(dirname "${annotation_fp}")
    else
        out_dir=$annotation_fp
    fi


    # Make transcriptome-index if directory doesn't exist
    if [ ! -d $out_dir\/$script_name-$script_ver ]; then
        mkdir $out_dir\/$script_name-$script_ver

        # For spike-ins first make merge genome file
        if [ -n "$spike_in" ]; then
            base_file=$(basename "${spike_in}")
            spike_in_fn=${base_file%.fa}
            bowtie2_index=$out_dir\/$script_name-$script_ver\/Bowtie2Index
            mkdir $bowtie2_index
            echo "* Generating merged genome, spike-in file ... "
            bowtie2-build -f $ucsc_genome_fn,$spike_in $bowtie2_index\/$ucsc_reference\_$spike_in_fn
            # Make sure the combined fa file is preserved so it dosen't get rebuilt
            bowtie2-inspect $bowtie2_index\/$ucsc_reference\_$spike_in_fn > $bowtie2_index\/$ucsc_reference\_$spike_in_fn.fa
        fi

        echo "* Reference file(s): '$ucsc_genome_fn'"
        echo "* Value of index_prefix: '$out_fn'"

        echo "* Generating transcriptome index... "
        tophat -G $annotation --transcriptome-index=$out_dir\/$script_name-$script_ver\/$out_fn $bowtie2_index

        # Get input and output files and then print out metadata.json file
        input_files=("$annotation" "$bowtie2_index")
        if [ -n "$spike_in" ]; then
            input_files+=("$spike_in")
        fi
        printf -v input "\"%s\"," "${input_files[@]}"
        input=${input%,}
        output_file=($out_dir\/$script_name-$script_ver/*)
        printf -v output "\"%s\"," "${output_file[@]}"
        output=${output%,}
        printf '{"script name":"%s","script version":"%s", "input files": [%s], "output files": [%s]}' "$script_name" "$script_ver" "$input" "$output" \
        | python -m json.tool > $out_dir\/$script_name-$script_ver\/$out_fn.metadata.json
        echo "* Finished."
    else
        echo "* Index has been generated"
    fi

}

main "$@"
