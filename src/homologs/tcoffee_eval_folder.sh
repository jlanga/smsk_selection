#!/usr/bin/env bash
set -euo pipefail

in_dir=$1
in_ext=$2
out_dir=$3
out_ext=$4
threads=$5


tcoffee_eval(){

    file_in=$1
    file_out=$2

    t_coffee \
        -infile "$file_in" \
        -evaluate \
        -output=score_ascii \
        -outfile "$file_out" \
        -n_core 1 \
        -quiet \
    || true

}

export -f tcoffee_eval

mkdir --parents "$out_dir"

find "$in_dir" -name "*.$in_ext" \
| sort -V \
| parallel --keep-order --jobs "$threads" \
    tcoffee_eval \
        "$in_dir/{/.}.$in_ext" \
        "$out_dir/{/.}.$out_ext" \
