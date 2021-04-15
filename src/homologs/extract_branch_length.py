#!/usr/bin/env python3

import argparse

from Bio import Phylo


def extract_branch_length(tree_fn, target_species):
    """Extract the dS/branch_length of the common ancestor of target species. 
    
    Missing species are not considered"""

    tree = Phylo.read(tree_fn, "newick")

    trx_to_species = {
        terminal.name: terminal.name.split("@")[0] for terminal in tree.get_terminals()
    }

    trx_in_tree = [  # Just names
        trx for trx, species in trx_to_species.items() if species in target_species
    ]

    leafs_in_tree = [  # Clades
        leaf for leaf in tree.get_terminals() if leaf.name in trx_in_tree
    ]

    branch_length = tree.common_ancestor(leafs_in_tree).branch_length

    print(branch_length)


def parse_arguments():
    """Argument parser for extract_branch_length"""
    parser = argparse.ArgumentParser(
        description="extract_branch_length.py: extract the branch length of the "
        " common ancestor of a set of species"
    )
    parser.add_argument(
        "-t",
        "--tree",
        help="Input species tree (Newick; not transcript tree)",
        required=True,
    )
    parser.add_argument(
        "-s",
        "--species",
        help="Target species, separated by commas. If some are missing, won't be taken into account",
        required=True,
    )
    return vars(parser.parse_args())


if __name__ == "__main__":

    ARGS = parse_arguments()

    extract_branch_length(ARGS["tree"], ARGS["species"].split(","))
