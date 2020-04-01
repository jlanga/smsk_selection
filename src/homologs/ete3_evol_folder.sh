#!/usr/bin/env bash

set -euo pipefail


ete3_evol(){

    tree=$1; shift
    msa=$1; shift
    image=$1; shift
    output=$1; shift
    species=$1; shift
    models=$*

    transcripts=$(python src/homologs/get_definers.py "$tree" "$species")

    ete3 evol \
        --resume \
        -t "$tree" \
        --alg "$msa" \
        --models ${models[@]} \
        --image "$image" \
        --mark "$transcripts" \
        --cpu 1 \
        --resume \
    > "$output" 

}

tree_folder=$1
msa_folder=$2
output_folder=$3
cores=$4
models=$5
species=$6

mkdir --parents "$output_folder"

export -f ete3_evol

find "$tree_folder" -name "*.nwk" \
| sort -V \
| parallel --keep-order --jobs "$cores" \
    ete3_evol \
        "$tree_folder/{/.}.nwk" \
        "$msa_folder/{/.}.fa" \
        "$output_folder/{/.}.pdf" \
        "$output_folder/{/.}.txt" \
        "$species" \
        "$models"

rm -rf /tmp/ete3-tmp