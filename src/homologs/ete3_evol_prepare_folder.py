#!/usr/bin/env python

"""ete3_evol_launcher.py
Script to prepare ete3 evol:
- slice and rename tree with sequences in alignment
- get the taxa1,,,taxa2 string that defines branches
- check that the alignment contains a minimum background and foreground taxa
- launch ete3 evol with all the required parameters
"""
from copy import deepcopy
import sys

from Bio import Phylo
from Bio import AlignIO

from helpers import \
    fix_dir_path, \
    process_folders

def get_species_to_seqid(alignment):
    """Get the species to seqid dictionary:
    
    {"aalo": "aalo@TRINITY_DN123937_c0_g1_i1.p1"}
    """
    return {
        sequence.name.split("@")[0]: sequence.name
        for sequence in alignment
    }


def get_seqid_to_species(alignment):
    """Get the seqid to species dictionary:
    
    {"aalo@TRINITY_DN123937_c0_g1_i1.p1": "aalo"}
    """
    return {
        sequence.name: sequence.name.split("@")[0]
        for sequence in alignment
    }


def keep_leafs(tree, leafs):
    """
    Trim the tree keeping only every leaf in `leafs`
    """
    tree_trimmed = deepcopy(tree)
    all_species = set(leaf.name for leaf in tree.get_terminals())
    for species_to_remove in all_species - set(leafs):
        tree_trimmed.prune(species_to_remove)
    return tree_trimmed


def rename_tree(tree, alignment):
    """Remove leafs not in alignment and rename them accordingly"""
    tree = deepcopy(tree)
    species_to_seqid = get_species_to_seqid(alignment)
    tree = keep_leafs(tree, species_to_seqid.keys())
    for leaf in tree.get_terminals():
        leaf.name = species_to_seqid[leaf.name]
    return tree


def has_enough_by_background_and_foreground(
        alignment, foreground_list, min_foreground=2, min_background=2
    ):
    """
    Return the alignment if it has at least min_foreground and min_background 
    sequences
    """
    alignment_ids = {
        sequence.id.split("@")[0]: sequence.id
        for sequence in alignment
    }
    n_foreground = len(set(foreground_list) & set(alignment_ids.keys()))
    n_background = len(set(alignment_ids.keys()) - set(foreground_list))
    if n_foreground >= min_foreground and n_background >= min_background:
        return True
    return False


def ete3_evol_prepare(
    tree_in_fn, alignment_in_fn,
    tree_out_fn,
    foreground_list, min_foreground=2, min_background=2
    ):
    """
    Read a species tree and alignment (nwk and fasta),
    Read the list of foreground taxa

    If there are enough foreground and background species:
        subset and rename the species tree into a protein tree
        write the alignment separately.
    """
    # Read inputs
    tree_in = Phylo.read(file=tree_in_fn, format="newick")
    alignment_in = AlignIO.read(handle=open(alignment_in_fn, "r"), format="fasta")

    # Slice and rename leafs in tree
    tree_out = rename_tree(tree=tree_in, alignment=alignment_in)

    # Check that there are enough sequences
    if has_enough_by_background_and_foreground(
            alignment_in, foreground_list, min_foreground, min_background
        ):
        Phylo.write(trees=tree_out, file=tree_out_fn, format="newick")


if __name__ == '__main__':

    if len(sys.argv) != 9:
        sys.stderr.write(
            "ERROR. Incorrect number of parameters. Example:\n"
            "ete3_evol_prepare.py species_tree.nwk msa_folder msa_extension "
            "output_folder output_extension human,chimp,bonobo 2 2\n"
        )
        sys.exit(1)

    TREE_IN_FN = sys.argv[1]
    MSA_DIR_IN = fix_dir_path(sys.argv[2])
    MSA_EXT_IN = sys.argv[3]
    TREE_DIR_OUT = fix_dir_path(sys.argv[4])
    TREE_EXT_OUT = sys.argv[5]
    FOREGROUND_LIST = sys.argv[6].split(",")
    MIN_FOREGROUND = int(sys.argv[7])
    MIN_BACKGROUND = int(sys.argv[8])

    def ete3_evol_prepare_wrapper(file_in, file_out):
        ete3_evol_prepare(
            tree_in_fn=TREE_IN_FN,
            alignment_in_fn=file_in,
            tree_out_fn=file_out,
            foreground_list=FOREGROUND_LIST,
            min_foreground=MIN_FOREGROUND,
            min_background=MIN_BACKGROUND
        )
    
    process_folders(
        MSA_DIR_IN, MSA_EXT_IN, TREE_DIR_OUT, TREE_EXT_OUT,
        ete3_evol_prepare_wrapper
    )



# if __name__ == '__main__':

#     if len(sys.argv) != 7:
#         sys.stderr.write(
#             "ERROR. Incorrect number of parameters. Example:\n"
#             "ete3_evol_prepare.py species_tree.nwk codon_alignment.fa "
#             "output.nwk human,chimp,bonobo 2 2\n"
#         )
#         sys.exit(1)

#     TREE_IN_FN = sys.argv[1]
#     ALIGNMENT_IN_FN = sys.argv[2]
#     TREE_OUT_FN = sys.argv[3]
#     FOREGROUND_LIST = sys.argv[4].split(",")
#     MIN_FOREGROUND = int(sys.argv[5])
#     MIN_BACKGROUND = int(sys.argv[6])

#     ete3_evol_prepare(
#         tree_in_fn=TREE_IN_FN,
#         alignment_in_fn=ALIGNMENT_IN_FN,
#         tree_out_fn=TREE_OUT_FN,
#         foreground_list=FOREGROUND_LIST,
#         min_foreground=MIN_FOREGROUND,
#         min_background=MIN_BACKGROUND
#     )
