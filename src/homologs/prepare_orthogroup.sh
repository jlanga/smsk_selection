#!/usr/bin/env bash

set -euo pipefail

prepare_orthogroup() {
    aln=$1
    cds=$2
    grep ^">" "$aln" | tr -d ">" | cut -f 1 -d " " | seqtk subseq "$cds" /dev/stdin
}



cds=$1
dir_in=${2%/}
ext_in=$3
dir_out=${4%/}
ext_out=$5

mkdir -p "$dir_out"

for alignment_file in "$dir_in"/*"{ext_in}"; do
    
    name=$(basename "${alignment_file}")
    cluster_id=${name%%.*}
    
    prepare_orthogroup \
        "${alignment_file}" \
        "$cds" \
    > "${dir_out}/${cluster_id}".aln.cds
    
    pal2nal.pl \
        "${dir_in}/${cluster_id}${ext_in}" \
        "${dir_out}/${cluster_id}".aln.cds \
        -output fasta \
        -nogap \
        -nomismatch \
    > "${dir_out}/${cluster_id}${ext_out}"
done
