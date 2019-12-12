#!/usr/bin/env python3

"""correct_tree_leaf_names.py: modify the sequence names from orthofinder
(species_transcript_id) to the names for the pdc3 scripts (species@transcript_id)
"""

import os
import sys
import re

from Bio import Phylo

def correct_tree_leaf_names(filename_in, filename_out):
    """Correct a single tree
    - Replaces the first _ into @: transition between orothofinder and pdc
    - Removes the ENA|id| since pdc removes it too
    """
    tree = Phylo.read(filename_in, "newick")
    ena_regex = re.compile(r'ENA\|[A-Z0-9]*\|')
    for terminal in tree.get_terminals():
        terminal.name = terminal.name.replace("_", "@", 1)
        terminal.name = ena_regex.sub("", terminal.name)
    Phylo.write(tree, filename_out, "newick")


if __name__ == '__main__':

    if len(sys.argv) != 5:
        sys.exit("correct_leaf_names.py in_dir .in_ext out_dir .out_ext")

    IN_DIR = sys.argv[1]
    IN_EXT = sys.argv[2]
    OUT_DIR = sys.argv[3]
    OUT_EXT = sys.argv[4]

    if not os.path.exists(OUT_DIR):
        os.makedirs(OUT_DIR)

    if not IN_DIR.endswith("/"):
        IN_DIR += "/"
    if not OUT_DIR.endswith("/"):
        OUT_DIR += "/"

    for in_file in os.listdir(IN_DIR):
        if in_file.endswith(IN_EXT):
            prefix = in_file.split("/")[-1].split(IN_EXT)[0]
            out_file = OUT_DIR + prefix + OUT_EXT
            correct_tree_leaf_names(IN_DIR + in_file, out_file)
