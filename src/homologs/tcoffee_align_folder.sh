#!/usr/bin/env bash
set -euo pipefail

in_dir=$1
in_ext=$2
out_dir=$3
out_ext=$4
threads=$5

tcoffee_align(){
    file_in=$1
    file_out=$2

    methods="muscle_msa mafftgins_msa t_coffee_msa kalign_msa"

    t_coffee \
        "$file_in" \
        -method "$methods" \
        -output=aln \
        -outfile "$file_out" \
        -n_core 1 \
        -quiet
}

export -f tcoffee_align

mkdir --parents "$out_dir"

find "$in_dir" -name "*.$in_ext" \
| sort -V \
| parallel --keep-order --jobs "$threads" \
    tcoffee_align \
        "$in_dir/{/.}.$in_ext" \
        "$out_dir/{/.}.$out_ext"
