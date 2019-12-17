#!/usr/bin/env bash
set -euo pipefail


in_dir=$1
in_ext=$2
out_dir=$3
out_ext=$4
threads=$5

mkdir -p "$out_dir"

for filein in "$in_dir"/*."$in_ext" ; do

    echo "$filein"

    cluster_id=$(basename -- "$filein" ."$in_ext")

    prefix="$out_dir"/"$cluster_id"
    fileout="$prefix"."$out_ext"

    if [[ -f "$fileout" ]] ; then

        >&2 echo "$fileout exits. Skipping."
    
    elif [[ $(grep -c ^">" "$filein") -ge 4 ]] ; then

        raxml-ng --check --msa "$filein" --model WAG

        raxml-ng \
            --msa "$filein" \
            --msa-format FASTA \
            --data-type AA \
            --prefix "$prefix" \
            --model WAG  \
            --seed 12345 \
            --threads "$threads" \
        1>&2

    else 

        >&2 echo "$filein has less than 4 sequences. Skipping."

    fi


done