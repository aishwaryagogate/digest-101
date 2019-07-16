#!/bin/bash
# merge-annotation.sh

script_name="merge-annotation.sh"
script_ver="1.1.0"

# Help function
usage() {
  echo "-h  --Help documentation for $script_name"
  echo "-r  --UCSC Reference genome (e.g. hg19, mm10)"
  echo "-a  --Additional Annotation file (e.g. LNCipedia.gtf)"
  echo "-s  --Spike-in set (e.g. ERCC)"
  echo "-o  --Path to the output directory"
  echo "-v  --Version of script"
  echo "Example: $script_name -r 'hg19' [-a '/path/to/gtf/lcnipedia.gtf /path/to/gtf/refseq.gtf'] [-s '/path/to/spike-in/ercc_92.fa.gz'] -o '/path/to/output'"
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
    module load bowtie2/gcc/2.2.3
    module load tophat/2.0.12
    module load iGenomes/2013-03-25
    module load python/2.7.x-anaconda
    module load cufflinks/2.2.1

    # Parsing options
    OPTIND=1 # Reset OPTIND
    while getopts :r:o:a:s:vh opt
        do
            case $opt in
                r) ucsc_reference=$OPTARG;;
                o) out_dir=$OPTARG;;
                a) additional_gtf=$OPTARG;;
                s) spike_in=$OPTARG;;
                v) version;;
                h) usage;;
            esac
        done

    shift $(($OPTIND -1))

    # Check for mandatory options
    if [[ -z $ucsc_reference ]] || [[ -z $out_dir ]]; then
        usage
    fi

    # Establish the GENCODE version to use and BOWTIEINDEX
    if [[ $ucsc_reference == 'hg19' ]]; then
        gencode_version='release_19'
        organism='human'
        annotation_fn='gencode.v19.annotation.gtf'
        ucsc_genome_fn=$iGENOMES_DB_DIR\/Homo_sapiens/UCSC/$ucsc_reference/Sequence/WholeGenomeFasta/genome.fa
    elif [[ $ucsc_reference == 'mm10' ]]; then
        gencode_version='release_M5'
        organism='mouse'
        annotation_fn='gencode.vM5.annotation.gtf'
        ucsc_genome_fn=$iGENOMES_DB_DIR\/Mus_musculus/UCSC/$ucsc_reference/Sequence/WholeGenomeFasta/genome.fa
    elif [[ $ucsc_reference == 'mm9' ]]; then
        gencode_version='release_M1'
        organism='mouse'
        annotation_fn='gencode.vM1.annotation.gtf'
        ucsc_genome_fn=$iGENOMES_DB_DIR\/Mus_musculus/UCSC/$ucsc_reference/Sequence/WholeGenomeFasta/genome.fa
    else
        usage
    fi

    # Make directory if doesn't exist
    gtf_dir=Gencode_$organism/$gencode_version
    if [ ! -e $out_dir\/$gtf_dir ]; then
        mkdir $out_dir\/$gtf_dir
    fi

    # Check if file exists, if not download file
    annotation_fp=$out_dir\/$gtf_dir\/$annotation_fn
    out_file_name=${annotation_fn%.gtf}
    if [ ! -e $annotation_fp ]; then
        echo "* Value of annotations: '$gencode_version'"
        echo "* Value of genome: '$ucsc_reference'"
        echo "* Value of organim: '$organism'"
        echo "* Downloading files..."
        wget ftp://ftp.sanger.ac.uk/pub/gencode/$gtf_dir\/$annotation_fn.gz -P $out_dir\/$gtf_dir\
        gunzip $out_dir\/$gtf_dir\/$annotation_fn.gz
    else
        echo "* Transciptome directory exists for '$organism','$ucsc_reference','$gencode_version'"
    fi

    # Make a single annotation file if additional files are provided
    if [ -n "$additional_gtf" ]; then
        for file in ${additional_gtf[@]}
        do
            additional_gtf_fn=$(basename "${file}")
            additional_gtf_prefix=${additional_gtf_fn%.gtf}
            gtf_dir+=\_$additional_gtf_prefix
            out_file_name+=\_$additional_gtf_prefix
        done

    fi

    # Add in spike-ins if provided
    if [ -n "$spike_in" ]; then
        gunzip $spike_in
        spike_in_pn=$(dirname "${spike_in}")
        spike_in_fn=$(basename "${spike_in}")
        spike_in_prefix=${spike_in_fn%.fa.gz}
        gtf_dir+=\_$spike_in_prefix
        out_file_name+=\_$spike_in_prefix
        spike_in_fp=$spike_in_pn\/$spike_in_prefix.fa
    fi

    # Make new directory if it doesn't exist
    if [ ! -d $out_dir\/$gtf_dir ]; then
        mkdir $out_dir\/$gtf_dir
    fi

    # Move files into directory
    if [ -n "$additional_gtf" ]; then
        for file in ${additional_gtf[@]}
        do
            mv $file $out_dir\/$gtf_dir
        done
    fi

    if [ -n "$spike_in" ]; then
        mv $spike_in_fp $out_dir\/$gtf_dir
    fi

    # Merge gtf files using awk script
    merge_path=$out_dir\/$gtf_dir
    cwd=`pwd`
    if [ ! -d $merge_path\/$script_name-$script_ver ]; then
        mkdir $merge_path\/$script_name-$script_ver

        # Get input files
        input_files=("$annotation_fp")
        if [ -n "$additional_gtf" ]; then
            input_files+=($merge_path/*gtf)
            printf '%s\n' "${input_files[@]}" > $merge_path\/GTF_list.txt
            echo "* Merge GTF files..."
            cd $merge_path
            cuffcompare -p 10 -s $ucsc_genome_fn -i GTF_list.txt
            cd $cwd
        fi

        if [ -n "$spike_in" ] && [ -n "$additional_gtf" ]; then
            input_files+=("$merge_path/${spike_in_fn%.gz}")
            inputs=($merge_path\/cuffcmp.combined.gtf,"$merge_path/${spike_in_fn%.gz}")
            echo "* Merge GTF files..."
            awk -f GTF.awk ${inputs[@]} > $merge_path\/$script_name-$script_ver\/$out_file_name.gtf
            mv $merge_path\/cuffcmp* $merge_path\/$script_name-$script_ver\/
        elif [ -n "$spike_in" ]; then
            input_files+=("$merge_path/${spike_in_fn%.gz}")
            echo "* Merge GTF files..."
            awk -f GTF.awk ${input_files[@]} > $merge_path\/$script_name-$script_ver\/$out_file_name.gtf
        else
            mv $merge_path\/cuffcmp.combined.gtf $merge_path\/$script_name-$script_ver\/$out_file_name.gtf
            mv $merge_path\/cuffcmp* $merge_path\/$script_name-$script_ver\/
        fi


        # Get input and output files and then print out metadata.json file
        printf -v input "\"%s\"," "${input_files[@]}"
        input=${input%,}
        output_file=($merge_path\/$script_name-$script_ver/*)
        printf -v output "\"%s\"," "${output_file[@]}"
        output=${output%,}
        printf '{"script name":"%s","script version":"%s", "input files": [%s], "output files": [%s]}' "$script_name" "$script_ver" "$input" "$output" | python -m json.tool > $merge_path\/$script_name-$script_ver/metadata.json
        echo "* Finished."
    else
        echo "* Index has been generated"
    fi

}

main "$@"
