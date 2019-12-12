#!/usr/bin/env python

"""extract_4d_sites_cds.py

Script to extract from a fasta file the 4-fold degenerate sites.

Relies on functions described in the Biopython cookbook:
https://biopython.org/wiki/Degenerated_Codons
"""

import sys

from Bio.Data.CodonTable import unambiguous_dna_by_id
from Bio import AlignIO


def altcodons(codon, table):
    """List codons that code for the same aminonacid / are also stop.

    @param codon
    @table code table id
    @return list of codons

    """
    tab = unambiguous_dna_by_id[table]

    if codon in tab.stop_codons:
        return tab.stop_codons

    try:
        aminoacid = tab.forward_table[codon]
    except KeyError:
        return []

    return [
        key for (key, value) in tab.forward_table.items()
        if value == aminoacid and key[0] == codon[0] and key[1] == codon[1]
    ]


def degeneration(codon, table):
    """Determine how many codons code for the same amino acid / are also stop

    @param codon the codon
    @param table code table id
    @param the number of codons also coding for the amino acid codon codes for
    """
    return len(altcodons(codon, table))


def is_x_degenerated(x_degeneracy, codon, table):
    """Determine if codon is x-fold degenerated.

    @param codon the codon
    @param table code table id
    @param true if x <= the degeneration of the codon
    """
    return x_degeneracy <= len(altcodons(codon, table))


def is_4d_or_gap(codon):
    """Check if the codon is 4d or a gap '---'"""
    if str(codon) == '---' or is_x_degenerated(4, codon, 1):
        return True
    return False



def extract_4d_sites(msa):
    """Subset the columns that are 4-fold degenerate sites

    @param msa the multiple sequence alignment
    @param the slice of the msa where all the columns are 4-fold degenerate sites
    """
    nt_length = msa.get_alignment_length()
    msa_degenerated = msa[:, :0]

    for codon_number in range(0, nt_length, 3):
        column_is_4d = all(
            is_4d_or_gap(str(record.seq[codon_number:codon_number+3]))
            for record in msa
        )
        if column_is_4d:
            msa_degenerated += msa[:, codon_number : codon_number + 3]
    return msa_degenerated



if __name__ == '__main__':

    if len(sys.argv) != 3:
        sys.stderr.write("ERROR: python3 extract_4d_sites_cds.py in.fa out.fa\n")
        sys.exit(1)

    FASTA_IN, FASTA_OUT = sys.argv[1], sys.argv[2]

    MSA_IN = AlignIO.read(handle=FASTA_IN, format="fasta")

    MSA_4D = extract_4d_sites(MSA_IN)
    if MSA_4D.get_alignment_length() > 0:
        AlignIO.write(alignments=MSA_4D, handle=FASTA_OUT, format="fasta")
    else:
        sys.stderr.write(FASTA_IN + ": No 4d sites were found\n")
