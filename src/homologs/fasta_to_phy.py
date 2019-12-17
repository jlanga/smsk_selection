#!/bin/env python

"""
fasta_to_phy
"""

import sys
import Bio.AlignIO


if __name__ == '__main__':

    if len(sys.argv) != 1:
        sys.exit(
            "ERROR! Incorrect number of parameters. Pipe a fasta to it and will "
            "print a .phy file:\n    fasta_to_phy.py < file.fasta > file.phy"
        )

    Bio.AlignIO.convert(
        sys.stdin, "fasta",
        sys.stdout, "phylip-sequential"
    )
