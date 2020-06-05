#!/usr/bin/env python3
"""shape_tree_from_msa.py: shape a species tree to the species and transcripts 
in the msa.
"""

import argparse

from Bio import \
    AlignIO, \
    Phylo


def shape_tree_from_msa(tree_in_fn, msa_in_fn, tree_out_fn):
    """Shape the tree as the msa:
    
    msa has sequences in the form species@transcript_id.

    - Remove from the tree each species not in the msa
    - Rename each species to the transcript id
    """

    tree = Phylo.read(tree_in_fn, "newick")
    msa = AlignIO.read(msa_in_fn, "fasta")

    # Get the spp - transcript dict
    species_to_trx = {
        name.split("@")[0]: name
        for name in [record.name for record in msa]
    }

    # Remove species
    species_in_tree = set(leaf.name for leaf in tree.get_terminals())
    species_not_in_msa = species_in_tree - set(species_to_trx.keys())
    for species in species_not_in_msa:
        tree.prune(species)

    # Rename the remaining
    for leaf in tree.get_terminals():
        leaf.name = species_to_trx[leaf.name]
    
    Phylo.write(trees=tree, file=tree_out_fn, format="newick")


def parse_arguments():
    """
    parser for run_fastcodeml.py
    """
    parser = argparse.ArgumentParser(
        description='shape_tree_as_msa.py: shape a species tree to the ones '
        'present in a msa file'
    )
    parser.add_argument(
        '-t', '--tree',
        help='Input species tree (Newick; not transcript tree)',
        required=True
    )
    parser.add_argument(
        '-m', '--msa',
        help='Input codon MSA (FASTA). Stops will be removed',
        required=True
    )
    parser.add_argument(
        '-o', '--output-tree',
        help='Output newick file',
        required=True
    )
    return vars(parser.parse_args())



if __name__ == '__main__':

    ARGS = parse_arguments()

    shape_tree_from_msa(
        tree_in_fn=ARGS["tree"],
        msa_in_fn=ARGS["msa"],
        tree_out_fn=ARGS["output_tree"]
    )
