#!/usr/bin/env bash
set -euo pipefail

# prepare_cds

parse_file(){
    file=$1
    species=$(basename "$file" | cut -f 1 -d ".")

    seqtk seq "$file" \
    | cut -f 1 -d " " \
    | paste - - \
    | sed "s/^>/>$species@/" \
    | sed 's/ENA|[A-Z0-9]*|//' \
    | tr "\t" "\n"
}

export -f parse_file

parallel --keep-order \
    parse_file {} \
::: "$@"