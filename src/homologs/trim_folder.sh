#!/usr/bin/env bash

in_dir=$1
in_ext=$2
out_dir=$3
out_ext=$4
threads=$5
min_occupancy=$6

mkdir --parents "$out_dir"

find "$in_dir" -name "*.$in_ext" \
| sort -V \
| parallel --keep-order --jobs "$threads" \
    pxclsq \
        --seqf "$in_dir/{/.}.$in_ext" \
        --outf "$out_dir/{/.}.$out_ext" \
        --prop "$min_occupancy"