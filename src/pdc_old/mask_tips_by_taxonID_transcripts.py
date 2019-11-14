#!/usr/bin/env python2

"""
Input is a dir of trees that end with ".tt". Change INTREE_FILE_ENDING if needed 

Mask both mono- and paraphyletic tips that belong to the same taxon
If only mask monophyletic tips, comment out this line:
curroot = paraphyly_masking(curroot,unamb_chrDICT)
Keep the tip that has the most un-ambiguous, well-aligned charactors in the trimmed alignment
"""

import newick3,phylo3,os,sys
from tree_utils import get_name,remove_kink
from seq import read_fasta_file

INTREE_FILE_ENDING = ".tt"

def get_clusterID(filename):
    """given a file name return the cluster id"""
    return filename.split(".")[0]
    
def mask_monophyletic_tips(curroot,unamb_chrDICT):
    going = True
    while going and len(curroot.leaves()) >= 4:
        going = False
        for node in curroot.iternodes(): #walk through nodes
            if not node.istip: continue #only look at tips
            for sister in node.get_sisters():
                if sister.istip and get_name(node.label)==get_name(sister.label): #masking
                    #print node.label,unamb_chrDICT[node.label],sister.label,unamb_chrDICT[sister.label]
                    if unamb_chrDICT[node.label] > unamb_chrDICT[sister.label]:
                        node = sister.prune()
                    else: node = node.prune()
                    if len(curroot.leaves()) >= 4:
                        if (node==curroot and node.nchildren==2) or (node!=curroot and node.nchildren==1):
                            node,curroot = remove_kink(node,curroot)
                    going = True
                    break
    return curroot
    
def mask_paraphyletic_tips(curroot,unamb_chrDICT):
    going = True
    while going and len(curroot.leaves()) >= 4:
        going = False
        for node in curroot.iternodes(): #walk through nodes
            if not node.istip: continue #only look at tips
            parent = node.parent
            if node == curroot or parent == curroot:
                continue #no paraphyletic tips for the root
            for para in parent.get_sisters():
                if para.istip and get_name(node.label)==get_name(para.label):
                    if unamb_chrDICT[node.label] > unamb_chrDICT[para.label]:
                        node = para.prune()
                    else: node = node.prune()
                    if len(curroot.leaves()) >= 4:
                        if (node==curroot and node.nchildren==2) or (node!=curroot and node.nchildren==1):
                            node,curroot = remove_kink(node,curroot)
                    going = True
                    break
    return curroot
    

def main(filename_in, filename_out, alignment_in, mask_para):
	with open(filename_in, "r") as filein, open(filename_out, "w") as fileout:
		intree = newick3.parse(filein.readline())
		clusterID = get_clusterID(filename_in)
		chrDICT = {}
		for sequence in read_fasta_file(alignment_in):
			for character in '-Xx?*':
				sequence.seq = sequence.seq.replace(character, "")  # ignore gaps, xs and Xs
			chrDICT[sequence.name] = len(sequence.seq)
		current_root = mask_monophyletic_tips(intree, chrDICT)
		if mask_para:
			currenst_root = mask_paraphyletic_tips(current_root, chrDICT)
		fileout.write(newick3.to_string(current_root)+"\n")


if __name__ == '__main__':
	
	if len(sys.argv) != 5:
		print(
			"python2 mask_tips_by_taxonID_transcripts.py in_tree out_tree "
			"alignment_in mask_paraphyletic(y/n)"
		)
		sys.exit(1)
	
	main(
		filename_in=sys.argv[1],
		filename_out=sys.argv[2],
		alignment_in=sys.argv[3],
		mask_para=sys.argv[4]
	)
