#!/usr/bin/env bash

in_dir=$1
in_ext=$2
out_dir=$3
out_ext=$4
threads=$5
min_occupancy=$6

mkdir --parents "$out_dir"

trim(){
    file_in=$1
    file_out=$2
    min_occupancy=$3

    pxclsq \
        --seqf "$file_in" \
        --outf /dev/stdout \
        --prop "$min_occupancy" \
    | seqtk seq \
    | paste - - \
    | awk '$2 != "-"' \
    | tr "\t" "\n" \
    > "$file_out"

    if [[ -s "$file_out" ]] ; then
        rm "$file_out"
    fi

}

find "$in_dir" -name "*.$in_ext" \
| sort -V \
| parallel --keep-order --jobs "$threads" \
    pxclsq \
        --seqf "$in_dir/{/.}.$in_ext" \
        --outf "$out_dir/{/.}.$out_ext" \
        --prop "$min_occupancy"