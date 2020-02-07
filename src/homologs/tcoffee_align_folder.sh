#!/usr/bin/env bash

set -euo pipefail

in_dir=$1
in_ext=$2
out_dir=$3
out_ext=$4
threads=$5


tcoffee_align(){
    pep=$1
    aln=$2

    methods="muscle_msa mafftgins_msa t_coffee_msa kalign_msa"

    if [[ ! -e "$aln" ]] ; then

        t_coffee \
            -seq "$pep" \
            -method "$methods" \
            -output=aln \
            -outfile "$aln" \
            -n_core 1 \
            -quiet || true
    
    fi
}

export -f tcoffee_align

mkdir --parents "$out_dir"

find "$in_dir" -name "*.$in_ext" -type f \
| sort --version-sort \
| parallel --jobs "$threads" \
    tcoffee_align \
        "$in_dir/{/.}.$in_ext" \
        "$out_dir/{/.}.$out_ext"

find . -maxdepth 1 -name "*.dnd" -delete