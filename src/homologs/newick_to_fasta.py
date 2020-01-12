#!/usr/bin/env python2

"""newick_to_fasta.py
Write to files the pep and cds fasta files contained in a newick tree
"""

import sys

from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio import Phylo


if __name__ == '__main__':

    if len(sys.argv) != 6:
        sys.stderr.write(
            "python2 newick_to_fasta.py tree.nwk in.cds in.pep out.cds out.pep")
        sys.exit(-1)

    # Parse arguments
    TREE_FN = sys.argv[1]
    CDS_IN_FN = sys.argv[2]
    PEP_IN_FN = sys.argv[3]
    CDS_OUT_FN = sys.argv[4]
    PEP_OUT_FN = sys.argv[5]

    # Parse tree
    try:
        TREE = Phylo.read(file=TREE_FN, format="newick")
        with open(CDS_IN_FN, "r") as cds_in, open(PEP_IN_FN, "r") as pep_in:
            CDS_DICT = {
                identifier.split()[0]: sequence
                for identifier, sequence in SimpleFastaParser(handle=cds_in)
            }
            PEP_DICT = {
                identifier.split()[0]: sequence
                for identifier, sequence in SimpleFastaParser(handle=pep_in)
            }

        # Write files
        with open(CDS_OUT_FN, "w") as cds_out, open(PEP_OUT_FN, "w") as pep_out:
            for terminal in TREE.get_terminals():
                identifier = terminal.name
                cds_out.write(">{identifier}\n{sequence}\n".format(
                    identifier=identifier,
                    sequence=CDS_DICT[identifier]
                ))
                pep_out.write(">{identifier}\n{sequence}\n".format(
                    identifier=identifier,
                    sequence=PEP_DICT[identifier]
                ))

    except ValueError:
        with open(CDS_OUT_FN, "w") as cds_out, open(PEP_OUT_FN, "w") as pep_out:
            cds_out.write("")
            pep_out.write("")
