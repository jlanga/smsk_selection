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
    in_ext=$2
    out_dir=$3
    
    in_dir=$(dirname "$file_in")
    cluster_id=$(basename -- "$file_in" ."$in_ext")

    prefix_out="$out_dir/$cluster_id"
    file_out="$out_dir/$cluster_id".raxml.bestTree


    if [[ -f "$file_out" ]] ; then

        >&2 echo "$file_out exits. Skipping."
    
    elif [[ $(grep -c ^">" "$file_in") -ge 4 ]] ; then

        raxml-ng --check --msa "$file_in" --model WAG

        raxml-ng \
            --msa "$file_in" \
            --msa-format FASTA \
            --data-type AA \
            --prefix "$prefix_out" \
            --model WAG  \
            --seed 12345 \
            --threads 1 \
        1>&2

    else 

        >&2 echo "$file_in has less than 4 sequences. Skipping."

    fi

}


export -f raxml_wrapper

parallel \
    --jobs "$threads" \
    raxml_wrapper {} "$in_ext" "$out_dir" \
::: "$in_dir"/*."$in_ext"

# for filein in "$in_dir"/*."$in_ext" ; do

#     echo "$filein"

#     cluster_id=$(basename -- "$filein" ."$in_ext")

#     prefix="$out_dir"/"$cluster_id"
#     fileout="$prefix"."$out_ext"

#     if [[ -f "$fileout" ]] ; then

#         >&2 echo "$fileout exits. Skipping."
    
#     elif [[ $(grep -c ^">" "$filein") -ge 4 ]] ; then

#         raxml-ng --check --msa "$filein" --model WAG

#         raxml-ng \
#             --msa "$filein" \
#             --msa-format FASTA \
#             --data-type AA \
#             --prefix "$prefix" \
#             --model WAG  \
#             --seed 12345 \
#             --threads "$threads" \
#         1>&2

#     else 

#         >&2 echo "$filein has less than 4 sequences. Skipping."

#     fi


# done