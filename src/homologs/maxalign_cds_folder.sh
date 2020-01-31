#!/usr/bin/env bash
set -euo pipefail

filter_dir=$1
filter_ext=$2
subset_dir=$3
subset_ext=$4
out_dir=$5
out_ext=$6
threads=$7


get_cds(){
    pep_in=$1
    cds_in=$2
    cds_out=$3

    perl src/maxalign.pl -p "$pep_in" "$cds_in" | seqtk seq - > "$cds_out" 

}

export -f get_cds

mkdir --parents "$out_dir"

find "$subset_dir" -name "*.$filter_ext" \
| sort --version-sort \
| parallel --keep-order --jobs "$threads" \
    get_cds \
        "$filter_dir/{/.}.$filter_ext" \
        "$subset_dir/{/.}.$subset_ext" \
        "$out_dir/{/.}.$out_ext"
