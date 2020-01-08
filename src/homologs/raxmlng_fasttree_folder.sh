#!/usr/bin/env bash
set -euo pipefail

#  bash src/homologs/raxmlng_folder.sh in_dir in_ext out_dir out_ext threads

in_dir=$1
in_ext=$2
out_dir=$3
threads=$4

mkdir -p "$out_dir"


raxml_wrapper() {
    
    file_in=$1
    file_out=$2

    in_dir=$(dirname "$file_in")
    in_ext=${file_in##*.}

    cluster_id=$(basename -- "$file_in" ."$in_ext")

    out_dir=$(dirname "$file_out")
  
    prefix_out="$out_dir/$cluster_id"
    file_out_raxml="$out_dir/$cluster_id".raxml.bestTree

    if [[ -f "$file_out" ]] ; then

        >&2 echo "$file_out exits. Skipping."
    
    elif [[ $(grep -c ^">" "$file_in") -ge 4 ]] ; then

        raxml-ng \
            --msa "$file_in" \
            --msa-format FASTA \
            --data-type AA \
            --prefix "$prefix_out" \
            --model WAG  \
            --seed 12345 \
            --threads 1 \
        1>&2

        ln "$file_out_raxml" "$file_out"
        rm "$prefix_out".raxml.{bestModel,log,mlTrees,rba,startTree}

    else 

        >&2 echo "$file_in has less than 4 sequences. Skipping."

    fi

}


fasttree_wrapper() {

    file_in=$1
    file_out=$2

    if [[ -f "$file_out" ]] ; then

        >&2 echo "$file_out exists. Skipping."

    else

        fasttree \
            -quiet \
            -wag \
            -out "$file_out" \
            "$file_in"
    
    fi 
}


conditional_execution () {

    file_in=$1
    file_out=$2

    n_seqs=$(grep -c ^">" "$file_in")

    if [[ "$n_seqs" -le 200 ]] ; then

        raxml_wrapper "$file_in" "$file_out"
    
    else 

        fasttree_wrapper "$file_in" "$file_out"

    fi

}




export -f \
    raxml_wrapper \
    fasttree_wrapper \
    conditional_execution

parallel \
    --jobs "$threads" \
    conditional_execution \
        "$in_dir/{/.}.$in_ext" \
        "$out_dir/{/.}.nwk" \
::: "$in_dir"/*."$in_ext"
