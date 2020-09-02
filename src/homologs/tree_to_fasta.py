#!/usr/bin/env python3

"""tree_to_fasta: convert a bunch of trees to fasta sequences, provided a fasta
file with all the sequences
"""

import sys
import os

from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio import Phylo


def fix_dir_path(path):
    """Add a / at the end if necessary"""
    return path if path.endswith("/") else path + "/"


def get_leafs(filein):
    """Just get the leaf names of a tree"""
    return [
        leaf_name.name
        for leaf_name in Phylo.read(filein, format="newick").get_terminals()
    ]


def fasta_to_dict(filein):
    """Parse fasta to a dict name: seq"""
    with open(filein, "r") as handle:
        return {
            identifier.split()[0]: sequence
            for identifier, sequence in SimpleFastaParser(handle)
        }


def process(fasta_dict, in_dir, in_ext, out_dir, out_ext):
    """Compose fasta clusters given the files in in_dir"""
    in_files = [
        in_dir + in_file for in_file in os.listdir(in_dir) if in_file.endswith(in_ext)
    ]

    for in_file in in_files:
        prefix = in_file.split("/")[-1].split(in_ext)[0]
        out_file = out_dir + prefix + out_ext
        sys.stderr.write(in_file + " => " + out_file + "\n")
        seq_names = get_leafs(in_file)
        with open(out_file, "w") as fout:
            for seq_name in seq_names:
                fout.write(
                    ">{seq_name}\n{seq}\n".format(
                        seq_name=seq_name, seq=fasta_dict[seq_name]
                    )
                )


if __name__ == "__main__":

    if len(sys.argv) != 6:
        sys.stderr.write(
            "ERROR: python2 tree_to_fasta.py all.fa in_dir in_ext out_dir out_ext"
        )

    FASTA = sys.argv[1]
    IN_DIR = fix_dir_path(sys.argv[2])
    IN_EXT = sys.argv[3]
    OUT_DIR = fix_dir_path(sys.argv[4])
    OUT_EXT = sys.argv[5]

    if not os.path.exists(OUT_DIR):
        os.makedirs(OUT_DIR)

    FASTA_DICT = fasta_to_dict(FASTA)

    process(FASTA_DICT, IN_DIR, IN_EXT, OUT_DIR, OUT_EXT)
