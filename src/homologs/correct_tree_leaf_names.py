#!/usr/bin/env python3

"""correct_tree_leaf_names.py: modify the sequence names from orthofinder
(species_transcript_id) to the names for the pdc3 scripts (species@transcript_id)
"""

import re
import sys

from Bio import Phylo

from helpers import fix_dir_path, process_folders


def correct_tree_leaf_names(filename_in, filename_out):
    """Correct a single tree
    - Replaces the first _ into @: transition between orothofinder and pdc
    - Removes the ENA|id| since pdc removes it too
    """
    tree = Phylo.read(filename_in, "newick")
    ena_regex = re.compile(r"ENA\|[A-Z0-9]*\|")
    for terminal in tree.get_terminals():
        terminal.name = terminal.name.replace("_", "@", 1)
        terminal.name = ena_regex.sub("", terminal.name)
    Phylo.write(tree, filename_out, "newick")


if __name__ == "__main__":

    if len(sys.argv) != 5:
        sys.exit("correct_leaf_names.py in_dir .in_ext out_dir .out_ext")

    IN_DIR = fix_dir_path(sys.argv[1])
    IN_EXT = sys.argv[2]
    OUT_DIR = fix_dir_path(sys.argv[3])
    OUT_EXT = sys.argv[4]

    process_folders(IN_DIR, IN_EXT, OUT_DIR, OUT_EXT, correct_tree_leaf_names)
