#!/bin/bash
# call-peaks-macs.sh

script_name="call-peaks-macs.sh"
script_ver="1.0.0"

main(){



    # Load required modules
    module load python/2.7.x-anaconda
    module load macs/1.4.2

    # Define sorted alignment file ending in .bed
    BEDFILE=$1

    # Define the input replicate 1 .bed
    INPUTBEDFILE1=$2

    #Define the input replicate 2 .bed
    INPUTBEDFILE2=$3

    BASEFILE=$(basename "${BEDFILE}")
    ONEPARENT=$(dirname "${BEDFILE}") # Get the first parent directory and go up
    OUTDIR=$(dirname "${ONEPARENT}")\/$script_name-$script_ver

    # Make Out directory if doesn't exist
    if [ ! -e $OUTDIR ]; then
        mkdir $OUTDIR
    fi

    # Combine Input replicates and sort if doesn't exist
    if [ ! -e $OUTDIR/Input_repCombined.sorted.BED ]; then

        # Combining replicates for peak calling:
        cat $INPUTBEDFILE1 $INPUTBEDFILE2 > $OUTDIR/Input_repCombined.BED

        # Sorting based on first two columns:
        sort -k1,1 -k2,2n $OUTDIR/Input_repCombined.BED > $OUTDIR/Input_repCombined.sorted.BED

        rm $OUTDIR/Input_repCombined.BED
    fi

    PRE=${BASEFILE%.BED}
    # Align if file doesn't exist
    if [ ! -e $OUTDIR/$PRE\_peaks.bed ]; then

        macs14 -t $BEDFILE -c $OUTDIR/Input_repCombined.sorted.BED -f BED -g hs -n $OUTDIR/$PRE

        # Get input and output files and then print out metadata.json file
        input_files=("$BEDFILE" "$INPUTBEDFILE1" "$INPUTBEDFILE2")
        printf -v input "\"%s\"," "${input_files[@]}"
        input=${input%,}
        output_file=($OUTDIR\/$PRE*)
        printf -v output "\"%s\"," "${output_file[@]}"
        output=${output%,}
        printf '{"script name":"%s","script version":"%s", "input files": [%s], "output files": [%s]}' "$script_name" "$script_ver" "$input" "$output"  | python -m json.tool > $OUTDIR/$PRE.metadata.json

    else
        echo "* $BASEFILE has peaks already called for it"
    fi
}
