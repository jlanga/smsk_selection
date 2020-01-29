#!/usr/bin/env bash
set -euo pipefail

filter_by_number_of_leafs() {

    file_in=$1
    file_out=$2
    min_leafs=$3

    n_leafs=$(python2.7 src/homologs/count_leafs.py < "$file_in")

    if [ "$n_leafs" -ge "$min_leafs" ] ; then
        ln "$file_in" "$file_out"
    fi

}

export -f filter_by_number_of_leafs

folder_in=$1
ext_in=$2
folder_out=$3
ext_out=$4
min_leafs=$5
threads=$6

mkdir --parents "$folder_out"

find "$folder_in" -name "*.$ext_in" \
| sort -V \
| parallel --jobs "$threads" --keep-order \
    filter_by_number_of_leafs \
        "$folder_in/{/.}.$ext_in" \
        "$folder_out/{/.}.$ext_out" \
        "$min_leafs" 
