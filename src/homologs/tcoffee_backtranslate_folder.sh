#!/usr/bin/env bash

set -euo pipefail

aln_dir=$1
aln_ext=$2
cds_dir=$3
cds_ext=$4
score_dir=$5
score_ext=$6
out_dir=$7
out_ext=$8
threads=$9

tcoffee_backtranslate(){
    aln=$1
    cds=$2
    score_ascii=$3
    fasta=$4

    if [[ ! -e "$fasta" ]] ; then

        HOME_4_TCOFFEE=$(pwd)
        export HOME_4_TCOFFEE

        t_coffee -other_pg seq_reformat \
            -in "$aln" \
            -in2 "$cds" \
            -struc_in "$score_ascii" \
            -struc_in_f number_aln \
            -output tcs_column_filter9_fasta \
            -out "$fasta" || true

    fi

}

export -f tcoffee_backtranslate

mkdir --parents "$out_dir"

find "$score_dir" -name "*.$score_ext" -type f \
| sort --version-sort \
| parallel --jobs "$threads" \
    tcoffee_backtranslate \
        "$aln_dir/{/.}.$aln_ext" \
        "$cds_dir/{/.}.$cds_ext" \
        "$score_dir/{/.}.$score_ext" \
        "$out_dir/{/.}.$out_ext"
