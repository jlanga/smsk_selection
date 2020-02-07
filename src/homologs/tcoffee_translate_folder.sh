#!/usr/bin/env bash

set -euo pipefail

in_dir=$1
in_ext=$2
out_dir=$3
out_ext=$4
threads=$5

tcoffee_translate(){
    cds=$1
    pep=$2

    if [[ ! -e "$pep" ]] ; then

        t_coffee -other_pg seq_reformat \
            -in "$cds" \
            -action +translate \
            -output fasta_seq \
            -out "$pep"
    
    fi

}


export -f tcoffee_translate

mkdir --parents "$out_dir"

find "$in_dir" -name "*.$in_ext" -type f \
| sort --version-sort \
| parallel --jobs "$threads" \
    tcoffee_translate \
        "$in_dir/{/.}.$in_ext" \
        "$out_dir/{/.}.$out_ext"

find . -maxdepth 1 -name "*.dnd" -delete