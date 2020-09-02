#!/usr/bin/env python

"""get_definers.py: script to get the two sequences that define the branch that 
define the elements specified
"""

import sys
from itertools import combinations

from Bio import Phylo

# def get_transcript_to_species(tree):
#     return {x.name: x.name.split("@")[0] for x in tree.get_terminals()}


def get_species_to_transcript(tree):
    """
    Get the {species: species@transcript_id} dictionary
    """
    return {x.name.split("@")[0]: x.name for x in tree.get_terminals()}


def translate_target_sequences(tree, target_species):
    """
    Translate the {species} set into {species@transcript_id} set
    """
    species_to_transcript = get_species_to_transcript(tree)
    return {
        species_to_transcript[species]
        for species in target_species
        if species in species_to_transcript
    }


def get_pair_that_defines_group(tree, target_transcripts):
    """Get a pair of elements in group that defines it as a subclade in tree
    
    tree must be a Phylo.tree object, and group a set of terminal nodes in the tree 
    
    return a list of two leafs or None if it not possible
    """
    for u, v in combinations(target_transcripts, 2):
        maximal_set = set(x.name for x in tree.common_ancestor(u, v).get_terminals())
        if target_transcripts.issubset(maximal_set):
            return [u, v]
    return None


def get_definers(tree_fn, target_species):
    """
    Get the two transcript_ids that define the specified branch in the tree from
    a set of target species.
    """
    tree = Phylo.read(tree_fn, "newick")

    target_transcripts = translate_target_sequences(tree, target_species)
    if len(target_transcripts) == 1:
        return target_transcripts

    return get_pair_that_defines_group(tree, target_transcripts)


if __name__ == "__main__":

    if len(sys.argv) != 3:
        sys.stderr.write("ERROR! python3 get_definers.py tree.nwk spp1,spp2,...,sppN\n")
        sys.exit(1)

    tree_fn = sys.argv[1]
    target_species = set(sys.argv[2].split(","))
    pair = get_definers(tree_fn, target_species)
    if pair is None:
        sys.stderr.write(f"All species are missing!: {target_species}\n")
        sys.stderr.write(f"Inspect tree in {tree_fn}\n")
        sys.exit(0)  # No error since it could break the pipeline
    if isinstance(pair, str):
        sys.stdout.write(pair)
    elif len(pair) != 2:
        sys.stderr.write(f"Not enough species: {pair}\n")
        sys.stderr.write(f"Target were: {target_species}\n")
        sys.stderr.write(f"Inspect tree in {tree_fn}\n")
        sys.stdout.write(list(pair)[0])
    else:
        sys.stdout.write(",,".join(pair) + "\n")
