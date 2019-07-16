#!/bin/bash
# call-transcripts-timecourse.sh

script_name="call-transcripts-timecourse.sh"
script_ver="1.0.0"

#Help function
usage() {
  echo "-h Help documentation for $script_name"
  echo "-f  --File path to Replicate 1 alignments"
  echo "-s  --File path to Replicate 2 alignments"
  echo "-r  --UCSC Reference genome (e.g. hg19, mm10)"
  echo "-o  --Path to output directory"
  echo "-e  --The experiment names."
  echo "-b  --The LtsProbB value"
  echo "-u  --The UTS value"
  echo "-v  --Version of script"
  echo "Example: $script_name -f 'foo_1.bam man_1.bam' -s 'foo_2.bam man_2.bam' -r 'hg19' -e 'foo man' -o '/path/to/output/dir/' -b -200 -u 5"
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
    module load R/3.1.0-intel

    # Parsing options
    OPTIND=1 # Reset OPTIND
    while getopts :f:s:r:e:b:u:o:vh opt
        do
            case $opt in
                f) aln1=$OPTARG;;
                s) aln2=$OPTARG;;
                r) ucsc_reference=$OPTARG;;
                e) exp=$OPTARG;;
                o) out=$OPTARG;;
                b) ltprob=$OPTARG;;
                u) uts=$OPTARG;;
                v) version;;
                h) usage;;
            esac
        done

    shift $(($OPTIND -1))

    # Check for mandatory options
    if [[ -z $aln1 ]] || [[ -z $aln2 ]] ||  [[ -z $ucsc_reference ]] || [[ -z $ltprob ]] || [[ -z $uts ]] || [[ -z $exp ]] || [[ -z $out ]]; then
        usage
    fi

    # Define the out directory
    out_dir=$out\/$script_name-$script_ver

    # Make sure directories exist
    if [ ! -e $out ]; then
        mkdir $out
    fi

    if [ ! -e $out_dir ]; then
        mkdir $out_dir
    fi

    # Run call-transcripts.R
    if [ ! -e $out_dir/metadata.json ]; then
        Rscript call-transcripts-timecourse.R --replicate1 $aln1 --replicate2 $aln2 --experiment $exp --genome $ucsc_reference --out $out_dir --ltprob $ltprob --uts $uts

        # Get input and output files and then print out metadata.json file
        input_files=("${array_aln1[@]}" "${array_aln2[@]}")
        printf -v input "\"%s\"," "${input_files[@]}"
        input=${input%,}
        output_file=($out_dir\/*)
        printf -v output "\"%s\"," "${output_file[@]}"
        output=${output%,}
        printf '{"script name":"%s","script version":"%s", "input files": [%s], "LtsProbB": [%s], "UTS": [%s], "output files": [%s]}' "$script_name" "$script_ver" "$input" "$ltprob" "$uts" "$output"  | python -m json.tool > $out_dir/metadata.json
    else
        aln1_fn=$(basename "$aln1")
        aln2_fn=$(basename "$aln2")
        echo "* HMM tune has been made "
    fi
}

main "$@"
