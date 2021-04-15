import sys,os
import phylo3,newick3
import tree_utils

"""
to change sequence names with taxa names and make trees more readable,
and output a new file named infile.names

Create a tabular file that each line contains
code	taxon_name
separated by tab
"""

if __name__ == "__main__":
	if len(sys.argv) != 3:
		print "python taxon_name_subst.py table treefile"
		sys.exit(0)
	
	DICT = {} #key is seq acronym, value is full taxon name, separated by tab
	with open(sys.argv[1], "rU") as infile:
		for line in infile:
			spls = line.strip().split("\t")
			if len(spls) > 1:
				DICT[spls[0]] = spls[1]
	print DICT
	
	#DIR = sys.argv[2]+"/"	
	#for i in os.listdir(DIR):
	treefile = sys.argv[2]
	"""
	#for alignments in fasta format
	infile = open(DIR+i,"r")
	outfile = open(DIR+i+".names","w")
	for line in infile:
		if line[0] == ">":
			if "@" in line: #for homolog alignments
				spls = (line[1:].strip()).split("@")
				taxonID, seqID = spls[0],spls[1]
				outfile.write('>'+DICT[taxonID]+"@"+seqID+"\n")
			else: #for ortho and species alignments
				id = line[1:].strip()
				if id in DICT:
					outfile.write('>'+DICT[id]+"\n")
				else: outfile.write(line)
		else: #no change
			outfile.write(line)
	infile.close()
	outfile.close()
	"""
	with open(treefile,"r") as infile:
		intree = newick3.parse(infile.readline())
	for i in intree.iternodes():
		if i.istip:
			print i.label
			if "@" in i.label: #for homolog trees
				spls = (i.label).split("@")
				taxonID, seqID = spls[0],spls[1]
				i.label = taxonID+"_"+DICT[taxonID]+"@"+seqID
			else: #for ortho and species trees with seqID removed
				try:
					i.label += "_"+DICT[i.label]
					#i.label = DICT[i.label]
				except:
					print i.lable,"not in the taxon table provided"

	with open(treefile+".name","w") as outfile:
		outfile.write(newick3.tostring(intree)+";\n")

