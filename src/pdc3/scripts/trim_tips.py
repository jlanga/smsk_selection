"""
Trim tips that sticking out (> relative_cutoff and >10 times longer than sister)
Also trim any tips that are > absolute_cutoff
"""

import os
import sys


import newick3
#from tree_utils import *
from tree_utils import remove_kink

#return the outlier tip, with abnormal high contrast and long branch
def check_countrast_outlier(node0, node1, above0, above1, relative_cutoff):
    """Check contrast outlier"""
    if node0.istip and above0 > relative_cutoff:
        if above1 == 0.0 or above0 / above1 > 10:
            return node0
    if node1.istip and above1 > relative_cutoff:
        if above0 == 0.0 or above1 / above0 > 10:
            return node1
    return None

def remove_a_tip(root, tip_node):
    """Remove a tip"""
    node = tip_node.prune()
    if len(root.leaves()) > 3:
        node, root = remove_kink(node, root)
        return root
    print("Less than four tips left")
    return None

def trim(curroot, relative_cutoff, absolute_cutoff):
    """Trim tree"""
    if curroot.nchildren == 2:
        _, _ = remove_kink(curroot, curroot)
    going = True
    while going and curroot is not None and len(curroot.leaves()) > 3:
        going = False
        for i in curroot.iternodes(order=1):  # POSTORDER
            if i.nchildren == 0:  # at the tip
                i.data['len'] = i.length
                if i.length > absolute_cutoff:
                    curroot = remove_a_tip(curroot, i)
                    going = True
                    break
            elif i.nchildren == 1:  # kink in tree
                remove_kink(i, curroot)
                going = True
                break
            elif i.nchildren == 2:  # normal bifurcating internal nodes
                child0, child1 = i.children[0], i.children[1]
                above0, above1 = child0.data['len'], child1.data['len']
                i.data['len'] = ((above0 + above1) / 2.) + i.length  #stepwise average
                outlier = check_countrast_outlier(child0, child1, above0, above1, relative_cutoff)
                if outlier is not None:
                    curroot = remove_a_tip(curroot, outlier)
                    going = True  #need to keep checking
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
                        child1, child2 = i.children[index1], i.children[index2]
                        above1, above2 = child1.data['len'], child2.data['len']
                        outlier = check_countrast_outlier(
                            child1, child2, above1, above2, relative_cutoff
                        )
                        if outlier is not None:
                            print(above1, above2)
                            curroot = remove_a_tip(curroot, outlier)
                            going = True  #need to keep checking
                            keep_checking = False  #to break the nested loop
                            break
                    if not keep_checking:
                        break
    return curroot

def main(folder, tree_file_ending, relative_cut, absolute_cut):
    """Process all trees in folder"""
    if folder[-1] != "/":
        folder += "/"
    filecount = 0
    for i in os.listdir(folder):
        if i.endswith(tree_file_ending):
            print(i)
            filecount += 1
            with open(folder + i, "r") as infile:
                intree = newick3.parse(infile.readline())
            outtree = trim(intree, float(relative_cut), float(absolute_cut))
            if outtree is not None:
                with open(folder + i +".tt", "w") as outfile:
                    outfile.write(newick3.tostring(outtree) + ";\n")
    assert filecount > 0, \
        "No file end with "+ tree_file_ending + " found in " + folder

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("python trim_tips.py DIR tree_file_ending relative_cutoff absolute_cutoff")
        sys.exit(0)

    FOLDER, TREE_FILE_ENDING, RELATIVE_CUT, ABSOLUTE_CUT = \
        sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]
    main(FOLDER, TREE_FILE_ENDING, RELATIVE_CUT, ABSOLUTE_CUT)
