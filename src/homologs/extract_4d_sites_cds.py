#!/usr/bin/env python

"""extract_4d_sites_cds.py

Script to extract from a fasta file the 4-fold degenerate sites.
"""

# pylint: disable=no-name-in-module

import sys

from Bio import AlignIO, Seq

from helpers import process_folders, fix_dir_path


def get_codon_dict():
    """
    Compose the codon to aminoacid dict: dict["AAA"] = "K"
    """
    codon_dict = {}
    bases = "ACGT"
    for base1 in bases:
        for base2 in bases:
            for base3 in bases:
                codon = base1 + base2 + base3
                codon_dict[codon] = Seq.translate(codon, table="Standard")
    return codon_dict


def get_fourfold_degenerate_list():
    """Get the list of fourfold degenerate codons"""
    fourfold_degenerate_codons = []
    codon_dict = get_codon_dict()
    for codon, aminoacid in codon_dict.items():
        first_two = codon[0:2]
        change = False
        for third_base in "ACGT":
            if Seq.translate(first_two + third_base) != aminoacid:
                change = True
        if change is False:
            fourfold_degenerate_codons.append(codon)
    return fourfold_degenerate_codons


# 4d sites are:
# ACN -> T,
# CCN -> P, CGN -> R, CTN -> L
# GCN -> A, GGN -> G, GTN -> V
# TCN -> S
# The first two letters indicate what to search for


def extract_fourfold_degenerate_sites(alignment):
    """
    Extract the fourfold degenerate sites from a CDS alignment.

    4d sites are the third base of a codon, of codons such that
    any change in that base produce the same aminoacid

    Such codons are:
    ACN -> T,
    CCN -> P, CGN -> R, CTN -> L
    GCN -> A, GGN -> G, GTN -> V
    TCN -> S

    It is enough to check the first two bases to know if a codon is degenerate
    """

    # Initialize data
    degenerated_codons = set(  # Let gaps pass
        ["AC", "CC", "CG", "CT", "GC", "GG", "GT", "TC", "--"]
    )

    n_columns = alignment.get_alignment_length()
    if n_columns % 3 != 0:
        sys.stderr.write(
            "ERROR: "
            "Number of columns in alignment is not multiple of three. "
            "I need codons."
        )
        sys.exit(1)
    n_codons = int(n_columns / 3)
    n_sequences = len(alignment)

    output_alignment = alignment[:, :0]  # Just sequence names

    for codon_index in range(n_codons):

        # Get column
        codon_column = [
            str(alignment[sequence_index][3 * codon_index : 3 * codon_index + 3].seq)
            for sequence_index in range(n_sequences)
        ]

        # Get the first two bases, since they determine if they are degenerate
        two_bases_set = set(codon[0:2] for codon in codon_column)

        # Check that all are degenerate or gap
        if not all(two_bases in degenerated_codons for two_bases in two_bases_set):
            continue

        # Check that they are the same two start bases or a gap
        if len(two_bases_set) == 1 or (
            len(two_bases_set) == 2 and "--" in two_bases_set
        ):
            output_alignment += alignment[:, 3 * codon_index + 2 : 3 * codon_index + 3]

    return output_alignment


def filter_4ds_sites(file_in, file_out):
    msa_in = AlignIO.read(handle=file_in, format="fasta")
    msa_4d = extract_fourfold_degenerate_sites(msa_in)
    if msa_4d.get_alignment_length() > 0:
        AlignIO.write(alignments=msa_4d, handle=file_out, format="fasta")
    else:
        sys.stderr.write(file_in + ": No fourfold-degenerate site was found\n")


if __name__ == "__main__":

    if len(sys.argv) != 5:
        sys.stderr.write(
            "ERROR: python3 extract_4d_sites_cds.py folder_in extension_in "
            "folder_out extension_out\n"
        )
        sys.exit(1)

    IN_DIR = fix_dir_path(sys.argv[1])
    IN_EXT = sys.argv[2]
    OUT_DIR = fix_dir_path(sys.argv[3])
    OUT_EXT = sys.argv[4]

    process_folders(IN_DIR, IN_EXT, OUT_DIR, OUT_EXT, filter_4ds_sites)
