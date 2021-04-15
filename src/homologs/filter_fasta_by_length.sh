#!/usr/bin/env bash

set -euo pipefail

# bash filter_fasta_by_length fasta_in length fasta_out

< "$1" \
paste - - \
| awk -v len="$2" 'length($2) >= len' \
| tr "\t" "\n" \
> "$3"
