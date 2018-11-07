#!/usr/bin/env python3

import sys
from Bio import SeqIO

if len(sys.argv) != 3:
    exit(
        """ERROR!: Incorrect number of arguments!
        Usage: python3 filter_longest_protein_per_gene.py fasta_out fasta_in
        - is allowed for /dev/stdin or /dev/stdout
        """
    )

fasta_out = sys.argv[1]
fasta_in = sys.argv[2]

if fasta_in == "-":
    fasta_in = sys.stdin

if fasta_out == "-":
    fasta_out = sys.stdout


def compute_gene_to_longest_protein(protein_dict):
    """Build a dict gene: protein, where protein is the longest protein that
    that gene has. In case of tie, the first seen is the one that remains.
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

    protein_dict = SeqIO.to_dict(SeqIO.parse(handle=fasta_in, format="fasta"))

    gene_to_protein = compute_gene_to_longest_protein(protein_dict)

    with open(fasta_out, "w") as f_out:
        SeqIO.write(
            sequences=gene_to_protein.values(),
            format="fasta",
            handle=f_out
        )
