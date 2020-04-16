#!/usr/bin/env bash

set -euo pipefail


ete3_evol(){

    tree=$1; shift
    msa=$1; shift
    # image=$1; shift
    codeml=$1; shift
    output=$1; shift
    species=$1; shift
    omega=$1; shift
    models=$*

    transcripts=$(python src/homologs/get_definers.py "$tree" "$species")

    #xvfb-run --auto-servernum \
    ete3 evol \
        --resume \
        -t "$tree" \
        --alg "$msa" \
        --models ${models[@]} \
        `# --image "$image"` \
        -o "$codeml" \
        --codeml_param "omega,$omega" \
        --mark "$transcripts" \
        --cpu 1 \
        --resume \
        --codeml_binary "$PWD/bin/codeml" \
    > "$output" 

}

tree_folder=$1
msa_folder=$2
output_folder=$3
cores=$4
models=$5
species=$6

mkdir --parents \
    `# "$output_folder/images"` \
    "$output_folder/values" \
    "$output_folder/codeml"

export -f ete3_evol

find $tree_folder -name '*.nwk' | sort -V \
| parallel -a - --keep-order --jobs "$cores" \
    ete3_evol \
        "$tree_folder/{1/.}.nwk" \
        "$msa_folder/{1/.}.fa" \
        `# "$output_folder/images/{1/.}.{2}.pdf"` \
        "$output_folder/codeml/" \
        "$output_folder/values/{1/.}.{2}.txt" \
        "$species" \
        "{2}" \
        "$models" \
::: 0.5 1.0 1.5
