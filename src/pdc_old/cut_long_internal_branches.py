"""
Cut long internal branches
If no long branch, copy the input tree to the out dir for the next round

The input trees can be .tre, .tre.mm or .tre.mm.tt depends on the workflow
"""
import sys

import newick3
from tree_utils import get_front_names, remove_kink


def count_taxa(node):
    """given a node count how many taxa it has in frount"""
    return len(set(get_front_names(node)))


def cut_long_internal_branches(curroot, cutoff):
    """cut long branches and output all subtrees with at least 4 tips"""
    going = True
    subtrees = [] #store all subtrees after cutting
    while going:
        going = False #only keep going if long branches were found during last round
        for node in curroot.iternodes(): #Walk through nodes
            if node.istip or node == curroot:
                continue
            if node.nchildren == 1:
                node, curroot = remove_kink(node, curroot)
                going = True
                break
            child0, child1 = node.children[0], node.children[1]
            if node.length > cutoff:
                print node.length
                if not child0.istip and \
                    not child1.istip and \
                    child0.length + child1.length > cutoff:
                    print child0.length + child1.length
                    if count_taxa(child0) >= 4:
                        subtrees.append(child0)
                    if count_taxa(child1) >= 4:
                        subtrees.append(child1)
                else:
                    subtrees.append(node)
                node = node.prune()
                if len(curroot.leaves()) > 2: #no kink if only two left
                    node, curroot = remove_kink(node, curroot)
                    going = True
                break
    if count_taxa(curroot) >= 4:
        subtrees.append(curroot) #write out the residue after cutting
    return subtrees


def main(
        masked_tree_fn, cutted_tree_fn, internal_branch_length_cutoff, min_taxa
    ):
    """cut long branches and output subtrees as .subtre files
    if uncut and nothing changed betwee .tre and .subtree
    copy the original .tre file to the outdir
    """
	# Read masked tree
    with open(masked_tree_fn, "r") as filein, \
        open(cutted_tree_fn, "w") as fileout:
        intree = newick3.parse(filein.readline())
        num_taxa = count_taxa(intree)
        min_taxa = int(min_taxa)
        if num_taxa < min_taxa:
            sys.stderr.write(
                "Tree has {num_taxa} less than {min_taxa} taxa\n".format(
                    num_taxa=num_taxa,
                    min_taxa=min_taxa
                )
            )
            exit(-1)
        else:
            subtrees = cut_long_internal_branches(intree, internal_branch_length_cutoff)
            if not subtrees:
                sys.stderr.write(
                    "No tree with at least " + min_taxa + " taxa"
                )
            else:
                count = 0
                outsizes = ""
                for subtree in subtrees:
                    if count_taxa(subtree) >= min_taxa:
                        if subtree.nchildren == 2:  # fix bifurcating roots from cutting
                            _, subtree = remove_kink(subtree, subtree)
                        count += 1
                        fileout.write(newick3.tostring(subtree) + "\n")
                        outsizes += str(len(subtree.leaves())) + ", "
                sys.stderr.write(
                    str(count) + " tree(s) written. Sizes: " + outsizes + "\n"
                )



if __name__ == "__main__":
    if len(sys.argv) not in [5, 6]:
        sys.stderr.write(
            "python cut_long_internal_branches.py intree_fn outtree_fn "
            "internal_branch_length_cutoff minimal_taxa\n"
        )
        sys.exit(1)

    main(
        masked_tree_fn=sys.argv[1],
        cutted_tree_fn=sys.argv[2],
        internal_branch_length_cutoff=sys.argv[3],
        min_taxa=sys.argv[4]
    )
