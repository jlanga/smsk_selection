#!/usr/bin/env python

"""
Filter a fasta alignment according to its occupancy:

filter_by_occupancy.py fasta_raw.fa fasta_trimmed.fa 0.5
"""

import os
import sys

from helpers import fasta_to_dict


def filter_by_occupancy(filename_in, filename_out, min_occupancy=0.5):
    """
    Filter an alignment in fasta format according to the occupancy of the 
    columns. Store the results in fasta format.
    """
    fasta_raw = fasta_to_dict(filename_in)
    
    n_sequences = len(fasta_raw.keys())
    alignment_length = len(fasta_raw[tuple(fasta_raw.keys())[0]])
    
    columns = tuple(
        "".join(
            fasta_raw[seqname][column_index] 
            for seqname in fasta_raw.keys()
        )
        for column_index in range(alignment_length)
    )
    
    columns_to_keep = []
    for column_number, column in enumerate(columns):
        n_gaps = column.count('-')
        if float(n_gaps) / float(n_sequences) <= 1 - min_occupancy:
            columns_to_keep.append(column_number)
    
    fasta_trimmed = {}
    for seqname, sequence in fasta_raw.items():
        fasta_trimmed[seqname] = "".join(
            fasta_raw[seqname][column_to_keep] 
            for column_to_keep in columns_to_keep
        )
    
    if not os.path.exists(os.path.dirname(filename_out)):
        os.makedirs(os.path.dirname(filename_out))

    with open(filename_out, "w") as f_out:
        for seqname, sequence in fasta_trimmed.items():
            f_out.write(
                ">{seqname}\n{sequence}\n".format(
                    seqname=seqname,
                    sequence=sequence
                )
            )



if __name__ == '__main__':

    if len(sys.argv) != 4:
        sys.stderr.write(
            "ERROR: incorrect number of arguments.\n"
            "python filter_by_occupancy.py fastain fastaout min_occupancy\n"
        )
        sys.exit(1)

    FASTA_IN = sys.argv[1]
    FASTA_OUT = sys.argv[2]
    MIN_OCCUPANCY = float(sys.argv[3])

    filter_by_occupancy(FASTA_IN, FASTA_OUT, MIN_OCCUPANCY)