#!/usr/bin/env python

"""
submodule to split a fasta file into multiple fasta files according to a
folder of fasta files.
Example: we have a folder with multiple protein files, and we want the cds of
every single msa as codons in a separate folder for post-processing.
"""

# pylint: disable=no-name-in-module

import sys

from helpers import \
    fasta_to_dict, \
    fix_dir_path, \
    process_folders


def subset_file(pep_fn, cds_fn, cds_dict):
    """Write to cds_fn the cds sequences that are in pep_fn"""
    with open(cds_fn, "w") as cds_out:
        for seqid in fasta_to_dict(pep_fn).keys():
            cds_out.write(
                ">{seqid}\n{sequence}\n".format(
                    seqid=seqid,
                    sequence=cds_dict[seqid]
                )
            )


if __name__ == '__main__':

    if len(sys.argv) != 6:
        sys.stderr.write(
            "Error. Usage: python split_cds.py folder_in ext_in folder_out "
            "ext_out all_cds.fa"
        )
        sys.exit(1)

    IN_DIR = fix_dir_path(sys.argv[1])
    IN_EXT = sys.argv[2]
    OUT_DIR = fix_dir_path(sys.argv[3])
    OUT_EXT = sys.argv[4]
    CDS_FN = sys.argv[5]

    CDS_DICT = fasta_to_dict(CDS_FN)


    def process(pep_fn, cds_fn):
        """Fix parameter"""
        return subset_file(pep_fn, cds_fn, CDS_DICT)

    process_folders(IN_DIR, IN_EXT, OUT_DIR, OUT_EXT, process)
    