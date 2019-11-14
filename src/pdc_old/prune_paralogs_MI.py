"""
Input: homolog trees
Output: individual orthologs trees

if a tip is longer than the LONG_TIP_CUTOFF
and also long than 10 times its sister, cut it off
This is to fix the leftover trees that frequently has some long tips in it

If not to output 1-to-1 orthologs, for example, already analysed these
set OUTPUT_1to1_ORTHOLOGS to False
"""
import newick3,phylo3,os,sys
import trim_tips
import tree_utils

OUTPUT_1to1_ORTHOLOGS = True 

def get_clusterID(filename):
	return filename.split(".")[0]
	
def get_front_score(node):
	front_labels = tree_utils.get_front_labels(node)
	num_labels = len(front_labels)
	num_taxa = len(set([tree_utils.get_name(i) for i in front_labels]))
	if num_taxa == num_labels:
		return num_taxa
	return -1
	
def get_back_score(node,root):
	back_labels = tree_utils.get_back_labels(node,root)
	num_labels = len(back_labels)
	num_taxa = len(set([tree_utils.get_name(i) for i in back_labels]))
	if num_taxa == num_labels:
		return num_taxa
	return -1
	
def prune(score_tuple,node,root,pp_trees):
	if score_tuple[0] > score_tuple[1]: #prune front
		print "prune front"
		pp_trees.append(node)
		par = node.prune()
		if par != None and len(root.leaves()) >= 3:
			par,root = tree_utils.remove_kink(par,root)
		return root,node == root
	else:
		if node != root: #prune back
			par = node.parent #par--node<
			par.remove_child(node)
			if par.parent != None:
				par,root = tree_utils.remove_kink(par,root)
		node.prune()
		print "prune back"
		pp_trees.append(root)
		if len(node.leaves()) >= 3:
			node,newroot = tree_utils.remove_kink(node,node)
		else:
			newroot = node
		return newroot,False #original root was cutoff, not done yet
			

if __name__ == "__main__":
	if len(sys.argv) != 7:
		print "python prune_paralogs_MI.py homoTreeDIR tree_file_ending relative_tip_cutoff absolute_tip_cutoff MIN_TAXA outDIR"
		print "LONG_TIP_CUTOFF is typically same value of the previous LONG_TIP_CUTOFF"
		sys.exit(0)

	inDIR = sys.argv[1]+"/"
	tree_file_ending = sys.argv[2]
	relative_tip_cutoff,absolute_tip_cutoff = float(sys.argv[3]),float(sys.argv[4])
	MIN_TAXA = int(sys.argv[5])
	outDIR = sys.argv[6]+"/"

	for i in os.listdir(inDIR):
		if not i.endswith(tree_file_ending): continue
		print i
		with open(inDIR+i,"r") as infile: #only 1 tree in each file
			intree = newick3.parse(infile.readline())
		curroot = intree
		pp_trees = []
		
		if get_front_score(curroot) >= MIN_TAXA: #No need to prune
			print "No pruning needed"
			if OUTPUT_1to1_ORTHOLOGS:
				os.system("cp "+inDIR+i+" "+outDIR+get_clusterID(i)+"_1to1ortho.tre")
		else: #scoring the tree
			going = True
			pp_trees = []
			
			while going: #python version of do..while loop
				highest = 0
				highest_node = None 
				score_hashes = {} #key is node, value is a tuple (front_score,back_score)
				for node in curroot.iternodes():
					front_score = get_front_score(node)
					back_score = get_back_score(node,curroot)
					score_hashes[node] = (front_score,back_score)
					if front_score > highest or back_score > highest:
						highest_node = node
						highest = max(front_score,back_score)
				if highest >= MIN_TAXA: #prune
					curroot,done = prune(score_hashes[highest_node],highest_node,curroot,pp_trees)
					if done or len(curroot.leaves()) < MIN_TAXA:
						going = False
						break
				else:
					going = False
					break
		
		if len(pp_trees) > 0:
			count = 1
			for tree in pp_trees:
				if tree.nchildren == 2:
					node,tree = trim_tips.remove_kink(tree,tree)
				tree = trim_tips.trim(tree,relative_tip_cutoff,absolute_tip_cutoff)
				if tree != None and len(tree.leaves()) >= MIN_TAXA:
					with open(outDIR+get_clusterID(i)+"_MIortho"+str(count)+".tre","w") as outfile:	
						outfile.write(newick3.tostring(tree)+";\n")
					count += 1
