#!/usr/bin/env bash
set -euo pipefail

in_dir=$1
in_ext=$2
out_dir=$3
out_ext=$4
threads=$5
cds_file=$6


subset(){
    in=$1
    out=$2
    cds=$3

    seqtk seq "$in" \
    | paste - - \
    | cut -f 1 \
    | tr -d ">" \
    | seqtk subseq "$cds" - \
    > "$out"

}

export -f subset

mkdir --parents "$out_dir"

find "$in_dir" -name "*.$in_ext" \
| sort --version-sort \
| parallel --keep-order --jobs "$threads" \
    subset \
        "$in_dir/{/.}.$in_ext" \
        "$out_dir/{/.}.$out_ext" \
        "$cds_file"