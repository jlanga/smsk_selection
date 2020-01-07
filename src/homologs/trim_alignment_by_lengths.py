#!/usr/bin/env python2

'''trim_alignment_by_length.py: go seq by seq seeing if the sequence has at
least min_length amino acids ("-" doesn't count)
'''

import sys

from Bio import SeqIO


if __name__ == '__main__':

    if len(sys.argv) != 3:
        sys.stderr.write(
            "ERROR: trim_alignment_by_lengths.py in.fa min_length"
        )

    FASTA_IN = sys.argv[1]
    MIN_LENGTH = int(sys.argv[2])

    with open(FASTA_IN, "r") as handle:
        for name, seq in SeqIO.FastaIO.SimpleFastaParser(handle):
            if len(seq.replace("-", "")) >= MIN_LENGTH:
                sys.stdout.write(
                    ">{name}\n{seq}\n".format(
                        name=name,
                        seq=seq
                    )
                )
