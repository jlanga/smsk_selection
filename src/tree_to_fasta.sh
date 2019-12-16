#!/usr/bin/env bash

set -euo pipefail

all_pep=$1
in_dir=$2
in_ext=$3
out_dir=$4
out_ext=$5

mkdir -p "$out_dir"

for file in "$in_dir"/*."$in_ext" ; do
    
    cluster_id=$(basename -- "$file" ."$in_ext")

    < "$in_dir"/"$cluster_id"."$in_ext" \
    python2 src/extract_leafs.py  \
    | seqtk subseq "$all_pep" - \
    > "$out_dir"/"$cluster_id"."$out_ext"

done