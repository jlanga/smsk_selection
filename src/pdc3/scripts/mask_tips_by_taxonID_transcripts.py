"""
Input is a dir of trees that end with ".tt". Change INTREE_FILE_ENDING if
needed

Mask both mono- and paraphyletic tips that belong to the same taxon
If only mask monophyletic tips, comment out this line:
curroot = paraphyly_masking(curroot,unamb_chr_dict)
Keep the tip that has the most un-ambiguous, well-aligned charactors in the
trimmed alignment
"""

import os
import sys

import newick3
from tree_utils import get_name, remove_kink
from seq import read_fasta_file

INTREE_FILE_ENDING = ".tt"


def get_cluster_id(filename):
    """given a file name return the cluster id"""
    return filename.split(".")[0]


def mask_monophyletic_tips(curroot, unamb_chr_dict):
    """Mask the monophyletic tips"""
    going = True
    while going and len(curroot.leaves()) >= 4:
        going = False
        for node in curroot.iternodes():  # walk through nodes
            if not node.istip:
                continue  # only look at tips
            for sister in node.get_sisters():
                if sister.istip and \
                    get_name(node.label) == get_name(sister.label): # masking
                    # print(node.label, unamb_chr_dict[node.label],
                    #       sister.label, unamb_chr_dict[sister.label])
                    if unamb_chr_dict[node.label] > unamb_chr_dict[sister.label]:
                        node = sister.prune()
                    else:
                        node = node.prune()
                    if len(curroot.leaves()) >= 4:
                        if (node == curroot and node.nchildren == 2) or \
                            (node != curroot and node.nchildren == 1):
                            node, curroot = remove_kink(node, curroot)
                    going = True
                    break
    return curroot


def mask_paraphyletic_tips(curroot, unamb_chr_dict):
    """Mask the paraphyletic tips"""
    going = True
    while going and len(curroot.leaves()) >= 4:
        going = False
        for node in curroot.iternodes():  # walk through nodes
            if not node.istip:
                continue  # only look at tips
            parent = node.parent
            if node == curroot or parent == curroot:
                continue  # no paraphyletic tips for the root
            for para in parent.get_sisters():
                if para.istip and get_name(node.label) == get_name(para.label):
                    if unamb_chr_dict[node.label] > unamb_chr_dict[para.label]:
                        node = para.prune()
                    else: node = node.prune()
                    if len(curroot.leaves()) >= 4:
                        if (node == curroot and node.nchildren == 2) or \
                            (node != curroot and node.nchildren == 1):
                            node, curroot = remove_kink(node, curroot)
                    going = True
                    break
    return curroot


def main(tree_dir, cln_dir, para, intree_file_ending=INTREE_FILE_ENDING):
    """Process all the trees in treeDir"""
    if tree_dir[-1] != "/":
        tree_dir += "/"
    if cln_dir[-1] != "/":
        cln_dir += "/"
    assert para in "yn", "mask paraphyletic tips? (y/n)"
    mask_para = (para == "y")
    filecount = 0

    filematch = {}  #key is cluster_id, value is the .aln-cln file
    
    # Get file names
    filenames = [
        filename for filename in os.listdir(cln_dir) 
        if filename.endswith(".aln-cln")
    ]

    for filename in filenames:
        cluster_id = get_cluster_id(filename)
        assert cluster_id not in filematch, \
            f"The cluster_id {cluster_id} repeats in {cln_dir}"
        filematch[cluster_id] = filename

    for i in os.listdir(tree_dir):
        if i.endswith(intree_file_ending):
            with open(tree_dir + i, "r") as infile:
                intree = newick3.parse(infile.readline())
            print(i)
            cluster_id = get_cluster_id(i)
            filecount += 1
            chr_dict = {} #key is seqid, value is number of unambiguous chrs
            for sequence in read_fasta_file(cln_dir + filematch[cluster_id]):
                for character in '-Xx?*':
                    sequence.seq = sequence.seq.replace(character, "") #ignore gaps, xs and Xs
                chr_dict[sequence.name] = len(sequence.seq)
            curroot = mask_monophyletic_tips(intree, chr_dict)
            if mask_para:
                curroot = mask_paraphyletic_tips(curroot, chr_dict)
            with open(tree_dir + i + ".mm", "w") as outfile:
                outfile.write(newick3.tostring(curroot) + ";\n")
    assert filecount > 0, \
        "No file ends with " + intree_file_ending + " found in " + tree_dir


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print(
            """
            python mask_tips_by_taxonID_transcripts.py tree_dir aln-cln_dir
            mask_paraphyletic(y/n)
            """
        )
        sys.exit(0)

    TREE_DIR = sys.argv[1]
    CLN_DIR = sys.argv[2]
    PARA = sys.argv[3]
    main(TREE_DIR, CLN_DIR, PARA)
