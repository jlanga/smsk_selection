#!/usr/bin/env python3

"""
filter_longest_protein_per_gene
"""

import sys
from Bio import SeqIO

if len(sys.argv) != 3:
    exit(
        """ERROR!: Incorrect number of arguments!
        Usage: python3 filter_longest_protein_per_gene.py FASTA_OUT FASTA_IN
        - is allowed for /dev/stdin or /dev/stdout
        """
    )

FASTA_OUT = sys.argv[1]
FASTA_IN = sys.argv[2]

if FASTA_IN == "-":
    FASTA_IN = sys.stdin

if FASTA_OUT == "-":
    FASTA_OUT = sys.stdout


def compute_gene_to_longest_protein(protein_dict):
    """Build a dict gene: protein, where protein is the longest protein that
    that gene has.

    In case of tie, the first seen is the one that remains.
    """
    gene_to_protein = {}

    for protein_id in protein_dict.keys():

        gene_id, _ = protein_dict[protein_id].description.split()[1].split("~~")

        if gene_id in gene_to_protein:
            len_protein_in_dict = len(protein_dict[gene_to_protein[gene_id].id])
            len_new_protein = len(protein_dict[protein_id])

            if len_new_protein > len_protein_in_dict:
                gene_to_protein[gene_id] = protein_dict[protein_id]

        else:
            gene_to_protein[gene_id] = protein_dict[protein_id]

    return gene_to_protein


if __name__ == '__main__':

    PROTEIN_DICT = SeqIO.to_dict(SeqIO.parse(handle=FASTA_IN, format="fasta"))

    GENE_TO_PROTEIN = compute_gene_to_longest_protein(PROTEIN_DICT)

    with open(FASTA_OUT, "w") as f_out:
        SeqIO.write(
            sequences=GENE_TO_PROTEIN.values(),
            format="fasta",
            handle=f_out
        )
