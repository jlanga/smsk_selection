#!/usr/bin/env python
"""extract_leafs.py

Extract the leafs from a newick file
"""

import sys

from Bio import Phylo

if __name__ == '__main__':
    if len(sys.argv) != 1:
        sys.stderr.write("ERROR: python extract_leafs.py < in.nwk")

    print(len(Phylo.read(sys.stdin, format='newick').get_terminals()))
