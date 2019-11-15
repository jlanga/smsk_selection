"""
Extract ingroup clades when outgroups present

Prepare a taxon code file with each line look like (separated by tab):
IN	taxonID1
IN	taxonID2
OUT	taxonID3
"""

import phylo3,newick3,os,sys
from tree_utils import *

if __name__ == "__main__":
	if len(sys.argv) != 7:
		print("python extract_clades.py inDIR treefileending outDIR MIN_INGROUP_TAXA TAXON_CODE phypartsDIR")
		sys.exit(0)
	
	inDIR = sys.argv[1]+"/"
	treefileending = sys.argv[2]
	outDIR = sys.argv[3]+"/"
	MIN_INGROUP_TAXA = int(sys.argv[4])
	taxon_code_file = sys.argv[5]
	phypartsDIR = sys.argv[6]

	INGROUPS = []
	OUTGROUPS = []
	with open(taxon_code_file,"r") as infile:
		for line in infile:
			if len(line) < 3: continue
			spls = line.strip().split("\t")
			if spls[0] == "IN":
				INGROUPS.append(spls[1])
			elif spls[0] == "OUT":
				OUTGROUPS.append(spls[1])
			else:
				print("Check TAXON_CODE file format")
				sys.exit()
	if len(set(INGROUPS) & set(OUTGROUPS)) > 0:
		print("Taxon ID",set(INGROUPS) & set(OUTGROUPS),"in both ingroups and outgroups")
		sys.exit(0)
	print(len(INGROUPS),"ingroup taxa and",len(OUTGROUPS),"outgroup taxa read")
	print("Ingroups:",INGROUPS)
	print("Outgroups:",OUTGROUPS)
	
	for treefile in os.listdir(inDIR):
		if not treefile.endswith(treefileending): continue
		with open(inDIR+treefile,"r") as infile:
			 intree = newick3.parse(infile.readline())
		curroot = intree
		all_names = get_front_names(curroot)
		print(treefile)
		
		#check taxonIDs
		ingroup_names = []
		outgroup_names = []
		for name in all_names:
			if name in INGROUPS:
				ingroup_names.append(name)
			elif name in OUTGROUPS:
				outgroup_names.append(name)
			else:
				print(name,"not in ingroups or outgroups")
				sys.exit()
		if len(set(ingroup_names)) < MIN_INGROUP_TAXA:
			print("not enough ingroup taxa in tree")
			continue
		if len(outgroup_names) == 0:
			print("No outgroup in tree")
			continue

		inclades = extract_rooted_ingroup_clades(curroot,INGROUPS,OUTGROUPS,MIN_INGROUP_TAXA)
		inclade_count = 0
		for inclade in inclades:
			inclade_count += 1
			with open(outDIR+treefile+"."+str(inclade_count),"w") as outfile:
				outfile.write(newick3.tostring(inclade)+";\n")
			for node in inclade.iternodes():
				if node.istip:
					node.label = get_name(node.label) # output multi-labeled tree for phyparts
			with open(phypartsDIR+treefile+"."+str(inclade_count),"w") as outfile:
				outfile.write(newick3.tostring(inclade)+";\n")
		print(inclade_count,"clades extracted")
