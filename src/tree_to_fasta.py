#!/usr/bin/env python3

"""tree_to_fasta: convert a bunch of trees to fasta sequences, provided a fasta
file with all the sequences
"""

import sys
import os

from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio import Phylo


def get_leafs(filein):
    """Just get the leaf names of a tree"""
    return [
        leaf_name
        for leaf_name in Phylo.read(filein, format='newick').get_terminals()
    ]

def fasta_to_dict(filein):
    """Parse fasta to a dict name: seq"""
    with open(filein, "r") as handle:
        return {
            identifier.split()[0]: sequence
            for identifier, sequence in SimpleFastaParser(handle)
        }

if __name__ == '__main__':

    if len(sys.argv) != 6:
        sys.stderr.write(
            "ERROR: python2 tree_to_fasta.py all.fa in_dir in_ext out_dir out_ext"
        )

    FASTA = sys.argv[1]
    IN_DIR = sys.argv[2]
    IN_EXT = sys.argv[3]
    OUT_DIR = sys.argv[4]
    OUT_EXT = sys.argv[5]

    if not os.path.exists(OUT_DIR):
        os.makedirs(OUT_DIR)

    if not IN_DIR.endswith("/"):
        IN_DIR += "/"
    if not OUT_DIR.endswith("/"):
        OUT_DIR += "/"

    FASTA_DICT = fasta_to_dict(FASTA)

    IN_FILES = [
        in_file
        for in_file in os.listdir(IN_DIR)
        if in_file.endswith(IN_EXT)
    ]

    for in_file in IN_FILES:
        prefix = in_file.split("/")[-1].split(IN_EXT)[0]
        out_file = OUT_DIR + prefix + OUT_EXT

        seq_names = get_leafs(in_file)
        with open(out_file, "w") as fout:
            for seq_name in seq_names:
                fout.write(
                    ">{seq_name}\n{seq}".format(
                        seq_name=seq_name,
                        seq=FASTA_DICT[seq_name]
                    )
                )
