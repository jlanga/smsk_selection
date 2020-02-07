#!/usr/bin/env bash

set -euo pipefail

in_dir=$1
in_ext=$2
out_dir=$3
out_ext=$4
threads=$5


tcoffee_evaluate(){
    aln=$1
    score_ascii=$2

    if [[ ! -e "$score_ascii" ]] ; then

        t_coffee -evaluate \
            -infile "$aln" \
            -output=score_ascii \
            -n_core 1 \
            -outfile "$score_ascii" \
            -quiet || true
    
    fi

}

export -f tcoffee_evaluate

mkdir --parents "$out_dir"

find "$in_dir" -name "*.$in_ext" -type f \
| sort --version-sort \
| parallel --jobs "$threads" \
    tcoffee_evaluate \
        "$in_dir/{/.}.$in_ext" \
        "$out_dir/{/.}.$out_ext"

find . -maxdepth 1 -name "*.dnd" -delete