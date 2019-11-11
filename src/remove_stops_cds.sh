#!/usr/bin/env bash

# remove_stops_cds.sh
# Remove stop codons from nucleotide fasta
# remove_stops_cds.sh < in.cds > out.cds


seqtk seq -C - \
| paste - - \
| sed -e "s/TAG$//g" -e "s/TGA$//g" -e "s/TAA$//g" \
| tr "\t" "\n"
