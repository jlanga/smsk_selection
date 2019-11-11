#!/usr/bin/env bash

# remove_stops_pep.sh
# Remove stop aa from protein fasta
# remove_stops_pep.sh < in.pep > out.pep

seqtk seq -C - \
| paste - - \
| sed -e "s/*$//g" \
| tr "\t" "\n"
