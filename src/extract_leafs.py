#!/usr/bin/env python
"""extract_leafs.py

Extract the leafs from a newick file
"""

import sys

from Bio import Phylo

if __name__ == '__main__':
    if len(sys.argv) != 1:
        sys.stderr.write("ERROR: python extract_leafs.py < in.nwk")

    for leaf_name in Phylo.read(sys.stdin, format='newick').get_terminals():
        sys.stdout.write(str(leaf_name) + "\n")
