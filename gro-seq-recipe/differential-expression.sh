#!/bin/bash
# differential-expression.sh

script_name="differential-expression.sh"
script_ver="1.0.0"

#Help function
usage() {
  echo "-h Help documentation for $script_name"\\n
  echo "-t1 --File path to Replicate 1 transcripts."
  echo "-t2 --File path to Replicate 2 transcripts."
  echo "-a1 --File path to Replicate 1 alignments."
  echo "-a2 --File path to Replicate 2 alignments."
  echo "-e  --The experiment names."
  echo "Example: $script_name -t1 'foo_rep1.bed man_rep1.bed' -t2 'foo_rep2.bed man_rep2.bed' -a1 'foo1.bam man1.bam' -a2 'foo2.bam man2.bam' -e 'foo man'"
  exit 1
}



main(){

    # Load required modules
    module load python/2.7.5
    module load R/3.1.0

    # Parsing options
    while getopts :t1:t2:a1:a2:e:h opt
        do
            case $opt in
                t1) TRX1=$OPTARG;;
                t2) TRX2=$OPTARG;;
                a1) ALN1=$OPTARG;;
                a2) ALN2=$OPTARG;;
                e) EXP=$OPTARG;;
                h) usage;;
            esac
        done

    shift $(($OPTIND -1))

    # Run make-signal.R
    Rscript differential-expression.R -t1 ${TRX1} -t2 ${TRX2} -a1 ${ALN1} -a2 ${ALN2} -e ${EXP}

}
