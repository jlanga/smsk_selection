#!/usr/bin/env bash
set -euo pipefail

in_dir=$1
in_ext=$2
out_dir=$3
out_ext=$4
threads=$5
min_length=$6


mkdir --parents "$out_dir"

find "$in_dir" -name "*.$in_ext" \
| sort --version-sort \
| parallel --keep-order --jobs "$threads" \
    bash src/homologs/filter_fasta_by_length.sh \
        "$in_dir/{/.}.$in_ext" \
        "$min_length" \
        "$out_dir/{/.}.$out_ext"
