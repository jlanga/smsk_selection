#!/usr/bin/env bash
set -euo pipefail

filter_dir=$1
filter_ext=$2
out_dir=$3
out_ext=$4
all_cds=$5
threads=$6


get_cds(){
    pep_in=$1
    cds_in=$2
    cds_out=$3

    pxaa2cdn -a "$pep_in" -n "$cds_in" | seqtk seq - > "$cds_out" 

}

export -f get_cds

mkdir --parents "$out_dir"

find "$filter_dir" -name "*.$filter_ext" \
| sort --version-sort \
| parallel --keep-order --jobs "$threads" \
    get_cds \
        "$filter_dir/{/.}.$filter_ext" \
        "$all_cds" \
        "$out_dir/{/.}.$out_ext"

# # Remove failed files
# (find "$out_dir" -name "*.$out_ext" -type f -print0 \
# | xargs --null grep --with-filename "fastafile" \
# | cut -f 1 -d ":" \
# | xargs --no-run-if-empty rm || true)