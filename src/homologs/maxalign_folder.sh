#!/usr/bin/env bash
set -euo pipefail

in_dir=$1
in_ext=$2
out_dir=$3
out_ext=$4
threads=$5


maxalign(){
    in=$1
    out=$2

    perl src/maxalign.pl -p "$in" | seqtk seq - > "$out"
}

export -f maxalign

mkdir --parents "$out_dir"

find "$in_dir" -name "*.$in_ext" \
| sort -V \
| parallel --keep-order --jobs "$threads" \
    maxalign \
        "$in_dir/{/.}.$in_ext" \
        "$out_dir/{/.}.$out_ext"
