#!/usr/bin/env bash

# remove_stops_cds.sh
# Remove stop codons from nucleotide fasta
# remove_stops_cds.sh < in.cds > out.cds


seqtk seq - \
| paste - - \
| sed -e "s/TAA$|TGA$T|AA$//g" \
| tr "\t" "\n"