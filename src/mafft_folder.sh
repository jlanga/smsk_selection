#!/usr/bin/env bash
set -euo pipefail


in_dir=$1
in_ext=$2
out_dir=$3
out_ext=$4
threads=$5


mkdir -p "$out_dir"

parallel \
    --jobs "$threads" \
    mafft --amino --genafpair --maxiterate 1000 {} ">" "$out_dir"/{/.}."$out_ext" \
::: "$in_dir"/*."$in_ext"

# Second round of dead processes, this time one by one multithreaded
find "$out_dir" -name "*.$out_ext" -size 0 \
| parallel \
    --pipe \
    --jobs 1 \
    mafft \
        --amino \
        --genafpair \
        --maxiterate 1000 \
        --thread $threads \
        "$in_dir"/{/.}."$in_ext" \
    ">" "$out_dir"/{/.}."$out_ext" 
