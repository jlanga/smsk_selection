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

tcoffee_backtranslate(){
    aln=$1
    cds=$2
    score_ascii=$3
    fasta=$4

    if [[ ! -e "$fasta" ]] ; then

        t_coffee -other_pg seq_reformat \
            -in "$aln" \
            -in2 "$cds" \
            -struc_in "$score_ascii" \
            -struc_in_f number_aln \
            -output tcs_column_filter9_fasta \
            -out "$fasta" || true

    fi

}


tcoffee(){
    cds=$1
    cds_trimmed=$2

    pep="$cds_trimmed".pep.tmp
    aln="$cds_trimmed".aln.tmp
    score_ascii="$cds_trimmed".score.tmp

    tcoffee_translate "$cds" "$pep"
    tcoffee_align "$pep" "$aln"
    tcoffee_evaluate "$aln" "$score_ascii"
    tcoffee_backtranslate "$aln" "$cds" "$score_ascii" "$cds_trimmed"

    rm --force "$pep" "$aln" "$score_ascii"

}


export -f tcoffee_translate
export -f tcoffee_align
export -f tcoffee_evaluate
export -f tcoffee_backtranslate
export -f tcoffee

mkdir --parents "$out_dir"

find "$in_dir" -name "*.$in_ext" -type f \
| sort --version-sort \
| parallel --jobs "$threads" \
    tcoffee \
        "$in_dir/{/.}.$in_ext" \
        "$out_dir/{/.}.$out_ext"

find . -maxdepth 1 -name "*.dnd" -delete