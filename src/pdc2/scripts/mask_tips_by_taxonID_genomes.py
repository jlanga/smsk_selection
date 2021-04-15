"""
Input is a dir of trees

Mask monophyletic tips that belong to the same taxon
Keep the tip with the shortest terminal branch lenth
"""

import newick3,phylo3,os,sys

#if taxon id pattern changes, change it here
def get_name(name):
	return name.split("@")[0]

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

def monophyly_masking_by_bl(curroot):
	going = True
	while going and len(curroot.leaves()) >= 4:
		going = False
		for node in curroot.iternodes(): #walk through nodes
			if not node.istip: continue #only look at tips
			for sister in node.get_sisters():
				if sister.istip and get_name(node.label)==get_name(sister.label):
					if node.length > sister.length:
						node = node.prune()			
					else: node = sister.prune()
					if len(curroot.leaves()) >= 4:
						if node.nchildren==1 or (node==curroot and node.nchildren==2):
							node,curroot = remove_kink(node,curroot)
							#no kink if the original node had more than 2 children
					going = True
					break
	return curroot

if __name__ == "__main__":
	if len(sys.argv) != 2:
		print "python mask_tips_by_taxonID_genomes.py DIR"
		sys.exit()

	DIR = sys.argv[1]+"/"
	filecount = 0
	for i in os.listdir(DIR):
		if i[-3:] == ".tt": #only mask trees that have tips trimmed
			with open(DIR+i,"r") as infile:
				intree = newick3.parse(infile.readline())
			print i
			filecount += 1
			with open(DIR+i+".mm","w") as outfile:
				outfile.write(newick3.tostring(monophyly_masking_by_bl(intree))+";\n")
				
	if filecount == 0:
		print "No file name with 'best' or 'tt' or 'fasttree' found in the treDIR"
