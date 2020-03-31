#!/usr/bin/env bash

set -euo pipefail

tree=$1
in_dir=$2
in_ext=$3
out_dir=$4
out_ext=$5
threads=$6
foreground=$7
min_foreground=$8
min_background=$9

mkdir --parents "$out_dir"


parallel \
    --jobs "$threads" \
    python3 src/homologs/ete3_evol_prepare.py \
        "$tree" \
        "$in_dir/{/.}.$in_ext" \
        "$out_dir/{/.}.$out_ext" \
        "$foreground" \
        "$min_foreground" \
        "$min_background" \
::: "$in_dir"/*."$in_ext"