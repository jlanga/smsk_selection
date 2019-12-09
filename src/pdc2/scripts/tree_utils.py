import phylo3,newick3
import os,sys

def get_name(label):
	"""Get taxonID from a tip label"""
	return label.split("@")[0]
	#return (label.split("@")[0]).split("_")[0]
	
def get_clusterID(filename):
	return filename.split(".")[0]

def get_front_labels(node):
	"""given a node, return a list of front tip labels"""
	leaves = node.leaves()
	return [i.label for i in leaves]

def get_back_labels(node,root):
	"""given a node, return a list of back tip labels"""
	all_labels = get_front_labels(root)
	front_labels = get_front_labels(node)
	return set(all_labels) - set(front_labels)
	
def get_front_names(node):
	"""given a node, return a list of front tip taxonIDs
	list may contain identical taxonIDs"""
	labels = get_front_labels(node)
	return [get_name(i) for i in labels]

def get_back_names(node,root):
	"""given a node, return a list of back tip taxonIDs
	list may contain identical taxonIDs"""
	back_labels = get_back_labels(node,root)
	return [get_name(i) for i in back_labels]


def remove_kink(node,curroot):
	"""
	smooth the kink created by prunning
	to prevent creating orphaned tips
	after prunning twice at the same node
	"""
	if node == curroot and curroot.nchildren == 2:
		#move the root away to an adjacent none-tip
		if curroot.children[0].istip: #the other child is not tip
			curroot = phylo3.reroot(curroot,curroot.children[1])
		else: curroot = phylo3.reroot(curroot,curroot.children[0])
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

def pass_boot_filter(node,min_ave_boot):
	"""check whether the average bootstrap value pass a cutoff"""
	total = 0.0
	count = 0
	for i in node.iternodes():
		if not i.istip and i.parent != None:
			total += float(i.label)
			count += 1
	if count == 0: #extracted clades with only two tips
		return True
	boot_average = total / float(count)
	print boot_average
	return boot_average >= float(min_ave_boot)

def get_ortho_from_rooted_inclade(inclade):
	"""
	input a rooted tree
	cut appart bifucating nodes when duplicated taxonIDs are detected
	"""
	assert inclade.nchildren == 2, "input clade not properly rooted"
	orthologs = [] #store ortho clades
	clades = [inclade]
	while True:
		newclades = [] #keep track of subclades generated in this round
		for clade in clades:
			num_taxa = len(set(get_front_names(clade)))
			num_tips = len((get_front_labels(clade)))
			if num_taxa == num_tips: #no taxon repeats
				orthologs.append(clade)
			else: #has duplicated taxa
				for node in clade.iternodes(order=0): #PREORDER, root to tip
					if node.istip: continue
					#traverse the tree from root to tip
					child0,child1 = node.children[0],node.children[1]
					name_set0 = set(get_front_names(child0))
					name_set1 = set(get_front_names(child1))
					if len(name_set0.intersection(name_set1)) > 0:
						if node == clade:
							newclades += [child0,child1] #break by bifid at the base
						elif len(name_set0) > len(name_set1): #cut the side with less taxa
							node.remove_child(child1)
							child1.prune()
							node,clade = remove_kink(node,clade) #no rerooting here
							newclades += [clade,child1]
						else:
							node.remove_child(child0)
							child0.prune()
							node,clade = remove_kink(node,clade) #no rerooting here
							newclades += [clade,child0]
						break
		if newclades == []: break
		clades = newclades
	return orthologs

def extract_rooted_ingroup_clades(root,ingroups,outgroups,min_ingroup_taxa):
	"""
	input a tree with ingroups and at least 1 outgroups
	output a list of rooted ingroup clades
	"""
	inclades = []
	while True:
		max_score,direction,max_node = 0,"",None
		for node in root.iternodes():
			front,back = 0,0
			front_names_set = set(get_front_names(node))
			for name in front_names_set:
				if name in outgroups:
					front = -1
					break
				elif name in ingroups: front += 1
				else: sys.exit("Check taxonID "+name)
			back_names_set = set(get_back_names(node,root))
			for name in back_names_set:
				if name in outgroups:
					back = -1
					break
				elif name in ingroups: back += 1
				else: sys.exit("Check taxonID "+name)
			if front > max_score:
				max_score,direction,max_node = front,"front",node
			if back > max_score:
				max_score,direction,max_node = back,"back",node
		#print max_score,direction
		if max_score >= min_ingroup_taxa:
			if direction == "front":
				inclades.append(max_node)
				kink = max_node.prune()
				if len(root.leaves()) > 3:
					newnode,root = remove_kink(kink,root)
				else: break
			elif direction == "back":
				par = max_node.parent
				par.remove_child(max_node)
				max_node.prune()
				inclades.append(phylo3.reroot(root,par))#flip dirction
				if len(max_node.leaves()) > 3:
					max_node,root = remove_kink(max_node,max_node)
				else: break
		else: break
	return inclades
