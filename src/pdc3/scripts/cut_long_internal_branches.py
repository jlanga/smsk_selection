"""
Cut long internal branches
If no long branch, copy the input tree to the out dir for the next round

The input trees can be .tre, .tre.mm or .tre.mm.tt depends on the workflow
"""

import newick3,phylo3,os,sys,math
from tree_utils import get_front_names,remove_kink,get_front_labels
from shutil import copy

def count_taxa(node):
	"""given a node count how many taxa it has in frount"""
	return len(set(get_front_names(node)))
	
def cut_long_internal_branches(curroot,cutoff):
	"""cut long branches and output all subtrees with at least 4 tips"""
	going = True
	subtrees = [] #store all subtrees after cutting
	while going:
		going = False #only keep going if long branches were found during last round
		for node in curroot.iternodes(): #Walk through nodes
			if node.istip or node == curroot: continue
			if node.nchildren == 1:
				node,curroot = remove_kink(node, curroot)
				going = True
				break
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
	if count_taxa(curroot) >= 4:
		subtrees.append(curroot) #write out the residue after cutting
	return subtrees

def main(inDIR,file_ending,branch_len_cutoff,min_taxa,outDIR):
	"""cut long branches and output subtrees as .subtre files
	if uncut and nothing changed betwee .tre and .subtree
	copy the original .tre file to the outdir"""
	if inDIR[-1] != "/": inDIR += "/"
	if outDIR[-1] != "/": outDIR += "/"
	min_taxa = int(min_taxa)
	filecount = 0
	cutoff = float(branch_len_cutoff)
	print("cutting branches longer than",cutoff)
	for i in os.listdir(inDIR):
		if not i.endswith(file_ending): continue
		print(i)
		filecount += 1
		with open(inDIR+i,"r") as infile: #only 1 tree in each file
			intree = newick3.parse(infile.readline())
		try:
			with open(inDIR+i[:i.find(".tre")]+".tre","r") as infile: #the original .tre
				raw_tree_size = len(get_front_labels(newick3.parse(infile.readline())))
		except: # did not refine this round. Use the .tre.tt.mm tree
			raw_tree_size = len(get_front_labels(intree))
		num_taxa = count_taxa(intree)
		if num_taxa < min_taxa:
			print("Tree has",num_taxa,"less than", min_taxa,"taxa")
		else:
			print(".tre:",raw_tree_size,"tips; "+file_ending+": "+str(len(get_front_labels(intree)))+" tips")
			subtrees = cut_long_internal_branches(intree,cutoff)
			if len(subtrees) == 0:
				print("No tree with at least", min_taxa, "taxa")
			#elif raw_tree_size == len(subtrees[0].leaves()):
				#copy(inDIR+i,outDIR+i)
				#print "written to out directory unchanged"
			else:
				count = 0
				outsizes = ""
				for subtree in subtrees:
					if count_taxa(subtree) >= min_taxa:
						if subtree.nchildren == 2: #fix bifurcating roots from cutting
							temp,subtree = remove_kink(subtree,subtree)
						count += 1
						with open(outDIR+i.split(".")[0]+"_"+str(count)+".subtree","w") as outfile:
							outfile.write(newick3.tostring(subtree)+";\n")
						outsizes += str(len(subtree.leaves()))+", "
				print(count,"tree(s) wirtten. Sizes:",outsizes)
	assert filecount > 0, "No file end with "+file_ending+" in "+inDIR
			
		
if __name__ == "__main__":
	if len(sys.argv) != 6:
		print("python cut_long_internal_branches.py inDIR tree_file_ending internal_branch_length_cutoff minimal_taxa outDIR")
		sys.exit(0)

	inDIR,file_ending,branch_len_cutoff,min_taxa,outDIR = sys.argv[1:]
	main(inDIR,file_ending,branch_len_cutoff,min_taxa,outDIR)

		
