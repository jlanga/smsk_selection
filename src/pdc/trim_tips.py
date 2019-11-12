#!/usr/bin/env python2
"""
Trim tips that sticking out (> relative_cutoff and >10 times longer than sister)
Also trim any tips that are > absolute_cutoff
"""

# import os
import sys
import newick3
# import phylo3

from tree_utils import *

#return the outlier tip, with abnormal high contrast and long branch
def check_countrast_outlier(node0,node1,above0,above1,relative_cutoff):
    if node0.istip and above0>relative_cutoff:
        if above1 == 0.0 or above0/above1 > 10:
            return node0
    if node1.istip and above1>relative_cutoff:
        if above0 == 0.0 or above1/above0 > 10 :
            return node1
    return None


def remove_a_tip(root,tip_node):
    node = tip_node.prune()
    if len(root.leaves()) > 3:
        node,root = remove_kink(node,root)
        return root
    else:
        print("Less than four tips left")
        return None


def trim(curroot,relative_cutoff,absolute_cutoff):
    if curroot.nchildren == 2:
        temp,root = remove_kink(curroot,curroot)
    going = True
    while going and curroot != None and len(curroot.leaves()) > 3:
        going = False
        for i in curroot.iternodes(order=1):# POSTORDER
            if i.nchildren == 0: # at the tip
                i.data['len'] = i.length
                if i.length > absolute_cutoff:
                    curroot = remove_a_tip(curroot,i)
                    going = True
                    break
            elif i.nchildren == 1: # kink in tree
                remove_kink(i,curroot)
                going = True
                break
            elif i.nchildren == 2: # normal bifurcating internal nodes
                child0,child1 = i.children[0],i.children[1]
                above0,above1 = child0.data['len'],child1.data['len']
                i.data['len'] = ((above0+above1)/2.)+i.length #stepwise average
                outlier = check_countrast_outlier(child0,child1,above0,above1,relative_cutoff)
                if outlier != None:
                    curroot = remove_a_tip(curroot,outlier)
                    going = True #need to keep checking
                    break
            else: #3 or more branches from this node. Pair-wise comparison
                total_len = 0
                nchild = i.nchildren
                for child in i.children:
                    total_len += child.data['len']
                i.data['len'] = total_len / float(i.nchildren)
                keep_checking = True
                for index1 in range(nchild): #do all the pairwise comparison
                    for index2 in range(nchild):
                        if index2 <= index1:
                            continue #avoid repeatedly checking a pair
                        child1,child2 = i.children[index1],i.children[index2]
                        above1,above2 = child1.data['len'], child2.data['len']
                        outlier = check_countrast_outlier(child1,child2,above1,above2,relative_cutoff)
                        if outlier != None:
                            print(above1, above2)
                            curroot = remove_a_tip(curroot,outlier)
                            going = True #need to keep checking
                            keep_checking = False #to break the nested loop
                            break
                    if not keep_checking: break
    return curroot


def main(in_tree, out_tree, relative_cutoff, absolute_cutoff):
    """Open file, apply trim trees, write results to another file
    """
    with open(in_tree, "r") as infile, open(out_tree, "w") as outfile:
        original_tree = newick3.parse(infile.readline())
        trimmed_tree = trim(
            curroot=original_tree,
            relative_cutoff=float(relative_cutoff),
            absolute_cutoff=float(absolute_cutoff)
        )
        outfile.write(newick3.tostring(trimmed_tree)+";\n")


if __name__ == "__main__":
    if len(sys.argv) != 5:
        print(
            "python trim_tips.py in_tree out_tree relative_cutoff "
            "absolute_cutoff"
        )
        sys.exit(1)

    main(
        in_tree=sys.argv[1],
        out_tree=sys.argv[2],
        relative_cutoff=sys.argv[3],
        absolute_cutoff=sys.argv[4]
    )
