#!/bin/bash
# align-bwa-se.sh

script_name="align-bwa-se.sh"
script_ver="1.1.0"

#Help function
usage() {
  echo "-h Help documentation for $script_name"
  echo "-f  --Path to fastq files a comma-separated list"
  echo "-r  --UCSC Reference genome (e.g. hg19, mm10)"
  echo "-o  --Path to output directory"
  echo "-v  --Version of script"
  echo "Example: $script_name -f 'foo_1.fastq.gz' -r 'hg19' [-o '/path/to/output/dir/']"
  echo "Example: $script_name -f 'foo_1.fastq.gz,fastq2.fastq.gz' -r 'hg19' [-o '/path/to/output/dir/']"
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
        index=$iGENOMES_DB_DIR\/Homo_sapiens/UCSC/$ucsc_reference/Sequence/BWAIndex/genome.fa
    elif [ $ucsc_reference = 'mm10' ]; then
        index=$iGENOMES_DB_DIR\/Mus_musculus/UCSC/$ucsc_reference/Sequence/BWAIndex/genome.fa
    elif [ $ucsc_reference = 'mm9' ]; then
        index=$iGENOMES_DB_DIR\/Mus_musculus/UCSC/$ucsc_reference/Sequence/BWAIndex/genome.fa
    elif [ $ucsc_reference = 'spike-in' ]; then
        index=/project/GCRB/shared/GRO-seq_spikein/SpikeIns.fa
    elif [ $ucsc_reference = 'rDNA' ]; then
        index=/project/GCRB/shared/GRO-seq_rDNA/u13369.fa
    else
        usage
    fi

    # Define the output directory, if none defined make the location relative to first file
    raw_file_array=(${raw_file//,/ })
    if [ -z $out ]; then
        out_dir=$(dirname "${raw_file_array[0]}")\/$script_name-$script_ver/
    else
        out_dir=$out\/$script_name-$script_ver
    fi

    if [ ! -d $out_dir ]; then
        mkdir $out_dir
    fi

    # Define the output file name, based on the first file
    raw_fn=$(basename "${raw_file_array[0]}")
    output_fp=${raw_fn%_*}

    # Align if file doesn't exist
    if [ ! -e $out_dir/$output_fp.sorted.bam ]; then

        # If array is more than 1 then make tmp directory to hold files to merge.
        if [[ ${#raw_file_array[@]} = 1 ]]; then
            echo "* Value of reads: '${raw_file_array[0]}'"
            echo "* Value of genome index : '$index'"

            echo "* Map reads..."
            bwa aln -t 8 $index ${raw_file_array[0]} > $out_dir/$output_fp.sai
            bwa samse -n 1 -f $out_dir/$output_fp.sam $index $out_dir/$output_fp.sai ${raw_file_array[0]}
            samtools view -b -S $out_dir/$output_fp.sam -o $out_dir/$output_fp.unsorted.bam
            samtools sort $out_dir/$output_fp.unsorted.bam $out_dir/$output_fp.sorted
            rm $out_dir/$output_fp.sam
            rm $out_dir/$output_fp.unsorted.bam
        else
            # Make temp directory
            mkdir $out_dir/tmp

            # For each file align the files
            for raw_file in ${raw_file_array[@]}
            do
                echo "* Value of reads: '$raw_file'"
                echo "* Value of genome index : '$index'"

                echo "* Map reads..."
                raw_fn=$(basename "$raw_file")
                bwa aln -t 8 $index $raw_file > $out_dir/tmp/$raw_fn.sai
                bwa samse -n 1 -f $out_dir/tmp/$raw_fn.sam $index $out_dir/tmp/$raw_fn.sai $raw_file
                samtools view -b -S $out_dir/tmp/$raw_fn.sam -o $out_dir/tmp/$raw_fn.unsorted.bam
                samtools sort $out_dir/tmp/$raw_fn.unsorted.bam $out_dir/tmp/$raw_fn.sorted
            done

            # Make new sam header
            sorted_bam_array=($out_dir/tmp/*sorted.bam)
            echo "* Making new header file..."
            samtools view -H $out_dir/tmp/$raw_fn.sorted.bam | grep -v "@PG" > $out_dir/tmp/new_header.sam
            samtools_version=samtools 2>&1 | grep "Version:" | cut -d' ' -f2
            new_PG="@PG\tID:samtools\tPN:samtools\tCL:"$samtools_version
            numb=1
            for sorted_bam in ${sorted_bam_array[@]}
            do
                samtools view -H tophat_out/accepted_hits.bam | grep "@PG" | \
                awk -v num=$numb 'BEGIN{OFS="\t"}{$2="ID:bwa_"num; print $0}' >> $out_dir/tmp/new_header.sam
                numb+=1
            done


            # Merge all aligned files and remove tmp directory
            samtools merge  $out_dir/$output_fp.sorted.bam $out_dir/tmp/*.sorted.bam
            rm -fr $out_dir/tmp
        fi

        # Get input and output files and then print out metadata.json file
        input_files=("$index" "${raw_file_array[@]}")
        printf -v input "\"%s\"," "${input_files[@]}"
        input=${input%,}
        output_file=($out_dir\/$output_fp*)
        printf -v output "\"%s\"," "${output_file[@]}"
        output=${output%,}
        printf '{"script name":"%s","script version":"%s", "input files": [%s], "output files": [%s]}' "$script_name" "$script_ver" "$input" "$output"  | python -m json.tool > $out_dir\/$output_fp.metadata.json

    else
        echo "* $output_fp has already been aligned"
    fi

}

main "$@"
