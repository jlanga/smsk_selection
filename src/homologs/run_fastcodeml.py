#!/usr/bin/python

"""run_fastcodeml.py: submodule to run fastcodeml with both null and alternative
hytopthesis, multiple omega_0 values, and return the tree, phy and tsv files 
with the p-values
"""

import argparse
from copy import deepcopy
import os
from subprocess import check_output
import sys

from Bio import Phylo
from Bio import AlignIO
from scipy.stats import chi2


def get_species_to_seqid(msa):
    """Get the species to seqid dictionary:
    
    {"aalo": "aalo@TRINITY_DN123937_c0_g1_i1.p1"}
    """
    return {
        sequence.name.split("@")[0]: sequence.name
        for sequence in msa
    }


def get_seqid_to_species(msa):
    """Get the seqid to species dictionary:
    
    {"aalo@TRINITY_DN123937_c0_g1_i1.p1": "aalo"}
    """
    return {
        sequence.name: sequence.name.split("@")[0]
        for sequence in msa
    }


def keep_leafs(tree, leafs):
    """
    Trim the tree keeping only every leaf in `leafs`
    """
    tree = deepcopy(tree)
    all_species = set(leaf.name for leaf in tree.get_terminals())
    for species_to_remove in all_species - set(leafs):
        tree.prune(species_to_remove)
    return tree


def tree_as_in_msa(tree, msa):
    """Remove leafs not in alignment and rename them accordingly"""
    tree = deepcopy(tree)
    species_to_seqid = get_species_to_seqid(msa)
    tree = keep_leafs(tree, species_to_seqid.keys())
    for leaf in tree.get_terminals():
        leaf.name = species_to_seqid[leaf.name]
    return tree


def msa_has_enough_by_background_and_foreground(
        msa, foreground_list, min_foreground=2, min_background=2
    ):
    """
    Return the alignment if it has at least min_foreground and min_background 
    sequences
    """
    ids = set(sequence.id.split("@")[0] for sequence in msa)
    foreground = set(foreground_list)
    n_foreground = len(foreground & ids)
    n_background = len(ids - foreground)
    if n_foreground >= min_foreground and n_background >= min_background:
        return True
    return False


def tree_has_enough_background_and_foreground(
    tree, foreground_list, min_foreground=2, min_background=2):
    """Return if the alignment has enoung background and foreground leafs"""
    leafs = set(leaf.name for leaf in tree.get_terminals())
    foreground = set(foreground_list)
    n_foreground = len(foreground & leafs)
    n_background = len(leafs - foreground)
    if n_background >= min_background and n_foreground >= min_foreground:
        return True
    return False


def clean_tree(tree):
    """Remove the confidence attribute in every branch of the tree"""
    tree = deepcopy(tree)
    for clade in tree.depths().keys():
        if hasattr(clade, 'confidence'):
            del clade.confidence
    return tree


def mark_branch(tree, target_leafs):
    """Mark in the tree the common ancestor of the species in target_leafs"""
    tree = deepcopy(tree)
    targets_in_tree = [
        target for target in target_leafs 
        if target in [leaf.name for leaf in tree.get_terminals()]
    ]
    if len(targets_in_tree) == 1:  # targets is only a leaf
        tree.common_ancestor(targets_in_tree).name += "#1"
    else:
        tree.common_ancestor(targets_in_tree).name = "#1"
    return tree


def remove_stops(msa):
    """Remove the stop codon in the alignment"""
    msa_trimmed = deepcopy(msa)
    for sequence in msa_trimmed:
        if str(sequence.seq[-3:]) in {"TAG", "TGA", "TAA"}:
            sequence.seq = sequence.seq[0:-3]
    return msa_trimmed


def clean_msa_names(msa):
    """Remove strange characters in sequence names (for fastcodeml)
    
    fastcodeml only accepts letters, numbers and some symbols
    """
    msa = deepcopy(msa)
    for sequence in msa:
        sequence.id = sequence.id.replace("@", "_", 1)
    return msa


def clean_tree_names(tree):
    """Remove strange characters in tree names (for fastcodeml)
    
    fastcodeml only accepts letters, numbers and some symbols
    """
    tree = deepcopy(tree)
    for leaf in tree.get_terminals():
        leaf.name = leaf.name.replace("@", "_", 1)
    return tree


def write_phy(msa, fileout):
    """Write the msa in a nonstadard phylip format: One line per record"""
    with open(fileout, "w") as f:
        nseq = len(msa)
        nt = len(msa[0].seq)
        f.write(f" {nseq} {nt}\n")
        for seq in msa:
            f.write(f"{seq.id} {seq.seq}\n")


def run_fastcodeml(
    input_tree_fn, input_msa_fn, out_prefix, hypothesis, omega_zero, 
    fastcodeml_binary="fast"):
    """Run fastcodeml and get the log likelihood"""
    command = [
        f"{fastcodeml_binary}",
        "--number-of-threads", "1", 
        "--only-hyp", str(hypothesis),
        "--init-param", f"w0={omega_zero}",
        "--init-param", f"p0=0.1",
        #"--branch-lengths-fixed",
        "--seed", "42",
        f"{out_prefix}.nwk", 
        f"{out_prefix}.phy"
    ]
    output = [
        line
        for line in check_output(command).decode().strip().split("\n")
        if line != ""
    ]
    lnl = float([line for line in output if line.startswith("LnL")][0].split()[1])
    return lnl


def write_pvalues(out_prefix, omegas={0.5, 1.0, 1.5}, fastcodeml_binary="fast"):
    """Run fastcodeml using multiple starting omegas and both hypothesis"""
    with open(f"{out_prefix}.tsv", "w") as tsv_out:
        
        tsv_out.write(f"orthogroup\tomega_zero\tlnL0\tlnL1\tl\tpvalue\n")
        
        for omega_zero in omegas:
            lnl0 = run_fastcodeml(
                input_tree_fn=f"{out_prefix}.nwk",
                input_msa_fn=f"{out_prefix}.phy",
                out_prefix=out_prefix,
                hypothesis=0,
                omega_zero=omega_zero,
                fastcodeml_binary=fastcodeml_binary
            )
            lnl1 = run_fastcodeml(
                input_tree_fn=f"{out_prefix}.nwk",
                input_msa_fn=f"{out_prefix}.phy",
                out_prefix=out_prefix,
                hypothesis=1,
                omega_zero=omega_zero,
                fastcodeml_binary=fastcodeml_binary
            )
            delta_l = 2 * abs(lnl1 -lnl0)
            pvalue = 1 - chi2.cdf(delta_l, 1)
            group = out_prefix.split("/")[-1]
            tsv_out.write(
                f"{group}\t{omega_zero}\t{lnl0}\t{lnl1}\t{delta_l}\t{pvalue}\n"
            )


def parse_arguments():
    """
    parser for run_fastcodeml.py
    """
    parser = argparse.ArgumentParser(
        description='run_fastcodeml.py: run fastcodeml with multiple omega_0 '
        'and get the p-values for each of them'
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
        '-w', '--omega-zeros',
        help='Starting omega values to start the ML optimization, separated by '
            'commas. Default: "0.5,1.0,1.5"',
        required=False,
        default='0.5,1.0,1.5'
    )
    parser.add_argument(
        '-s', '--target-species',
        help='Target species for the branch-site test, as in the input tree, '
            'separated by commas. Missing species will be ommited',
        required=True
    )
    parser.add_argument(
        '-b', '--min-background',
        help='Minimum number of species in the background of the tree. Default: 2',
        required=False,
        default=2,
        type=int
    )
    parser.add_argument(
        '-f', '--min-foreground',
        help='Minimum number of species in the foreground of the tree. Default: 2',
        required=False,
        default=2,
        type=int
    )
    parser.add_argument(
        '-o', '--output-prefix',
        help='Output prefix where to store the tree and msa for fastcodeml and '
            'the output tsv. Containing folders must exist',
        required=True
    )
    parser.add_argument(
        '-c', '--fastcodeml-binary',
        help='Path to the fastcodeml binary. Default: "fast"',
        required=False,
        default='fast',
    )

    return vars(parser.parse_args())


if __name__ == '__main__':

    ARGS = parse_arguments()
    # print(ARGS)
    MSA_FN = ARGS["msa"]

    TREE = Phylo.read(ARGS["tree"], "newick")
    MSA = AlignIO.read(MSA_FN, "fasta")
    TARGET_SPP = set(ARGS["target_species"].split(","))
    OUT_PREFIX = ARGS["output_prefix"]
    OMEGAS = set(ARGS["omega_zeros"].split(","))

    if not msa_has_enough_by_background_and_foreground(
            msa=MSA,
            foreground_list=TARGET_SPP,
            min_background=ARGS["min_background"],
            min_foreground=ARGS["min_foreground"]
        ):
        sys.stderr.write(
            f"{MSA_FN}: Not enough species background or foreground species\n"
        )
        sys.exit(0)

    TREE = tree_as_in_msa(TREE, MSA)
    TREE = clean_tree(TREE)
    SPP_TO_TRX = get_species_to_seqid(MSA)
    TARGET_TRX = [SPP_TO_TRX[SPP] for SPP in TARGET_SPP if SPP in SPP_TO_TRX]
    TREE = mark_branch(TREE, TARGET_TRX)    
    TREE = clean_tree_names(TREE)

    MSA = remove_stops(MSA)
    MSA = clean_msa_names(MSA)

    Phylo.write(TREE, f"{OUT_PREFIX}.nwk", "newick")
    write_phy(MSA, f"{OUT_PREFIX}.phy")

    write_pvalues(
        out_prefix=OUT_PREFIX,
        omegas=OMEGAS,
        fastcodeml_binary=ARGS["fastcodeml_binary"]
    )

