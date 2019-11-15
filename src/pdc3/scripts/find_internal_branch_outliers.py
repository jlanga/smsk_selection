"""
Separate subtrees connected by long branches
Input trees: tree files from tip masking

compare if the tree after cutting has exactly the same number of tips as the
starting tree of the last round
copy the alignment and tree to the outdir for the next round
"""

import newick3,phylo3,os,sys,math
from Bio import SeqIO

#if taxon id pattern changes, change it here
def get_name(name):
	return name.split("@")[0]
	
def get_leaf_labels(leaves):
	labels = []
	for i in leaves:
		labels.append(i.label)
	return labels

def count_taxa(node):
	labels = get_leaf_labels(node.leaves())
	names = []
	for label in labels:
		names.append(get_name(label))
	return len(set(names))
	
#smooth the kink created by prunning
#to prevent creating orphaned tips after prunning twice at the same node
def remove_kink(node,curroot):
	if node == curroot and curroot.nchildren == 2:
		#move the root away to an adjacent none-tip internal node
		if curroot.children[0].istip: #the other child is not tip
			curroot = phylo3.reroot(curroot,curroot.children[1])
		else: #tree has >=4 leaves so the other node cannot be tip
			curroot = phylo3.reroot(curroot,curroot.children[0])
	#---node---< all nodes should have one child only now
	length = node.length + (node.children[0]).length
	par = node.parent
	kink = node
	node = node.children[0]
	#parent--kink---node<
	par.remove_child(kink)
	par.add_child(node)
	node.length = length
	return node,curroot

#cut long branches and output all subtrees regardless of size
def cut_long_internal_branches(curroot,cutoff):
	going = True
	subtrees = [] #store all subtrees after cutting
	while going:
		going = False #only keep going if long branches were found during last round
		for node in curroot.iternodes(): #Walk through nodes
			if node.istip or node == curroot: continue
			child0,child1 = node.children[0],node.children[1]
			if node.length > cutoff:
				print(node.length)
				if not child0.istip and not child1.istip and child0.length+child1.length>cutoff:
					print(child0.length + child1.length)
					if count_taxa(child0) >= 4:
						subtrees.append(child0)
					if count_taxa(child1) >= 4:
						subtrees.append(child1)						
				else: subtrees.append(node)
				node = node.prune()
				if len(curroot.leaves()) > 2: #no kink if only two left
					node,curroot = remove_kink(node,curroot)
					going = True
				break
	if count_taxa(curroot) >= min_taxa:
		subtrees.append(curroot) #write out the residue after cutting
	return subtrees
	

if __name__ == "__main__":
	if len(sys.argv) != 5:
		print("python find_internal_branch_outliers.py inDIR file_ending minimal_taxa outDIR")
		sys.exit(0)

	inDIR = sys.argv[1]+"/"
	file_ending = sys.argv[2]
	min_taxa = int(sys.argv[3])
	outDIR = sys.argv[4]+"/"
	
	filecount = 0
	l = len(file_ending)
	for i in os.listdir(inDIR):
		if i[-l:] != file_ending: continue
		filecount += 1
		clusterID = i.split(".")[0]
		
		#input tree would either look like treeID.tt.mm or treeID.tt (No '.' in treeID)
		with open(inDIR+i,"r") as infile: #only 1 tree in each file
			intree = newick3.parse(infile.readline())
		curroot = intree
		
		internal_branch_lengths = []
		for node in curroot.iternodes():
			if node.istip or node == curroot: continue
			internal_branch_lengths.append(node.length)
		internal_branch_lengths.sort()
		with open(inDIR+i+".brlen","w") as outfile:
			for lenth in internal_branch_lengths:
				outfile.write(str(lenth)+"\n")
	
	if filecount == 0:
		print("No file end with",file_ending,"found")
			
		
		
