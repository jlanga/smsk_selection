#!/usr/bin/env bash
set -euo pipefail


in_dir=$1
in_ext=$2
out_dir=$3
out_ext=$4
threads=$5

mkdir -p "$out_dir"

raxmlng_reduce() {
    
    file_in=$1
    file_out=$2

    raxml-ng --check --msa "$file_in" --model WAG
    rm "$file_in".raxml.log

    reduced="$file_in".raxml.reduced.phy

    if [[ -s "$reduced" ]] ; then

        # phy to fa
        tail -n+2 "$reduced" | sed 's/^/>/g' | tr " " "\n" > "$file_out"
        rm -f "$reduced"

    else

        cp "$file_in" "$file_out"

    fi

}


export -f raxmlng_reduce

parallel \
    --jobs "$threads" \
    raxmlng_reduce \
        "$in_dir/{/.}.$in_ext" \
        "$out_dir/{/.}.$out_ext" \
::: "$in_dir"/*."$in_ext"
