#!/bin/bash
# align-bwa-pe.sh

script_name="align-bwa-pe.sh"
script_ver="1.1.0"

#Help function
usage() {
  echo "-h Help documentation for $script_name"
  echo "-f  --Path to read1 fastq files a comma-separated list"
  echo "-s  --Path to read2 fastq files a comma-separated list"
  echo "-r  --UCSC Reference genome (e.g. hg19, mm10)"
  echo "-o  --Path to output directory"
  echo "-v  --Version of script"
  echo "Example: $script_name -f 'foo_1.fastq.gz' -s 'foo_2.fastq.gz -r 'hg19' [-o '/path/to/output/dir/']"
  echo "Example: $script_name -f 'foo_1.fastq.gz,fastq2_1.fastq.gz' -s 'foo_2.fastq.gz,fastq2_2.fastq.gz' -r 'hg19' [-o '/path/to/output/dir/']"
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
    while getopts :f:s:r:o:vh opt
        do
            case $opt in
                f) raw_file_1=$OPTARG;;
                s) raw_file_2=$OPTARG;;
                r) ucsc_reference=$OPTARG;;
                o) out=$OPTARG;;
                v) version;;
                h) usage;;
            esac
        done

    shift $(($OPTIND -1))

    # Check for mandatory options
    if [[ -z $raw_file_1 ]] || [[ -z $raw_file_2 ]]  || [[ -z $ucsc_reference ]]; then
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
    raw_file_array_1=(${raw_file_1//,/ })
    raw_file_array_2=(${raw_file_2//,/ })
    if [ -z $out ]; then
        out_dir=$(dirname "${raw_file_array_1[0]}")\/$script_name-$script_ver/
    else
        out_dir=$out\/$script_name-$script_ver
    fi

    if [[(${#raw_file_array_1[@]} -ne ${#raw_file_array_2[@]}) ]]; then
        echo "The number of arguments are not equal."
        exit 1
    fi

    if [ ! -d $out_dir ]; then
        mkdir $out_dir
    fi

    # Define the output file name, based on the first file
    raw_fn=$(basename "${raw_file_array_1[0]}")
    output_fp=${raw_fn%_*}

    # Align if file doesn't exist
    if [ ! -e $out_dir/$output_fp.sorted.bam ]; then

        # If array is more than 1 then make tmp directory to hold files to merge.
        if [[ ${#raw_file_array_1[@]} = 1 ]]; then
            echo "* Value of reads: '${raw_file_array_1[0]}, ${raw_file_array_2[0]}'"
            echo "* Value of genome index : '$index'"

            echo "* Map reads..."
            bwa aln -t 8 $index ${raw_file_array_1[0]}  > $out_dir/$output_fp\_1.sai
            bwa aln -t 8 $index ${raw_file_array_2[0]}  > $out_dir/$output_fp\_2.sai
            bwa sampe -n 1 $index  $out_dir/$output_fp\_1.sai $out_dir/$output_fp\_2.sai ${raw_file_array_1[0]} ${raw_file_array_2[0]} > $out_dir/$output_fp.sam
            samtools view -b -S $out_dir/$output_fp.sam -o $out_dir/$output_fp.unsorted.bam
            samtools sort $out_dir/$output_fp.unsorted.bam $out_dir/$output_fp.sorted
            rm $out_dir/$output_fp.sam
            rm $out_dir/$output_fp.unsorted.bam
        else
            # Make temp directory
            mkdir $out_dir/tmp

            # For each file align the files
            for i in ${!raw_file_array_1[@]}
            do
                echo "* Value of reads: '${raw_file_array_1[$i]}', ${raw_file_array_2[$i]}"
                echo "* Value of genome index : '$index'"

                echo "* Map reads..."
                raw_fn=$(basename "${raw_file_array_1[$i]}")
                bwa aln -t 8 $index ${raw_file_array_1[$i]}  > $out_dir/$raw_fn\_1.sai
                bwa aln -t 8 $index ${raw_file_array_2[$i]} > $out_dir/$raw_fn\_2.sai
                bwa sampe -n 1 $index  $out_dir/$raw_fn\_1.sai $out_dir/$raw_fn\_2.sai ${raw_file_array_1[$i]} ${raw_file_array_2[$i]} $out_dir/$raw_fn.sam
                samtools view -b -S $out_dir/$output_fp.sam -o $out_dir/$output_fp.unsorted.bam
                samtools sort $out_dir/$raw_fn.unsorted.bam $out_dir/$raw_fn.sorted
                rm $out_dir/$raw_fn.sam
                rm $out_dir/$raw_fn.unsorted.bam
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
        input_files=("$index" "${raw_file_array_1[@]}" "${raw_file_array_2[@]}")
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
