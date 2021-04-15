#!/usr/bin/env bash
set -euo pipefail

dir_in=$1
ext_in=$2
dir_out=$3
ext_out=$4


remove_gaps(){
    fasta_in=$1
    fasta_out=$2

    seqtk seq "$fasta_in" \
    | paste - - \
    | awk '{gsub(/-/,"",$2); print $1"\t"$2}' \
    | tr "\t" "\n" \
    > "$fasta_out"

}

export -f remove_gaps

mkdir --parents "$dir_out"

find "$dir_in" -name "*.$ext_in" -type f \
| sort --version-sort \
| parallel --jobs 1 \
    remove_gaps \
        "$dir_in/{/.}.$ext_in" \
        "$dir_out/{/.}.$ext_out"