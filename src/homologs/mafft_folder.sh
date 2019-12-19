#!/usr/bin/env bash
set -euo pipefail


in_dir=$1
in_ext=$2
out_dir=$3
out_ext=$4
threads=$5


mkdir -p "$out_dir"

run_mafft_conditionally() {

    file=$1
    threads=$2

    nseq=$(grep -c ^">" "$file")
    naminoacids=$(
        seqtk seq "$file" \
        | paste - - \
        | awk "BEGIN{min=0} {if (length($2) > min) min=length($2)} END{print min}"
    )

    if [[ $nseq -lt 200 ]] && [[ $naminoacids -lt 2000 ]] ; then

        mafft \
            --amino \
            --genafpair \
            --maxiterate 1000 \
            --thread "$threads" \
            --quiet \
            "$file"

    else

        mafft \
            --amino \
            --auto \
            --thread "$threads" \
            --quiet \
            "$file"
    fi

}

export -f run_mafft_conditionally

parallel \
    --jobs "$threads" \
    run_mafft_conditionally {} 1 ">" "$out_dir/{/.}.$out_ext" \
::: "$in_dir"/*."$in_ext"
