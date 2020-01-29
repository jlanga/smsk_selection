#!/usr/bin/env bash

set -euo pipefail

tcoffee_filter(){
    aln=$1
    eval=$2
    out=$3

    python src/homologs/tcoffee_filter.py \
        "$aln" "$eval" "$out"

}

export -f tcoffee_filter

aln_dir=$1
aln_ext=$2
eval_dir=$3
eval_ext=$4
out_dir=$5
out_ext=$6
threads=$7


mkdir --parents "$out_dir"

find "$eval_dir" -name "*.$eval_ext" \
| sort --version-sort \
| parallel --keep-order --jobs "$threads" \
    tcoffee_filter \
        "$aln_dir/{/.}.$aln_ext" \
        "$eval_dir/{/.}.$eval_ext" \
        "$out_dir/{/.}.$out_ext"
