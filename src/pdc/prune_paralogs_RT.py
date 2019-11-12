"""
If no outgroup, only output ingroup clades with no taxon repeats
If outgroup present, extract rooted ingroup clades and prune paralogs

Prepare a taxon file, with each line look like (separated by tab):
IN	taxonID1
IN	taxonID2
OUT	taxonID3
"""

import phylo3,newick3,os,sys
import tree_utils

def RT(homoDIR,tree_file_eneding,outDIR,min_ingroup_taxa,taxon_code_file_file):
	if homoDIR[-1] != "/": homoDIR += "/"
	if outDIR[-1] != "/": outDIR += "/"
	min_ingroup_taxa = int(min_ingroup_taxa)

	INGROUPS = []
	OUTGROUPS = []
	with open(taxon_code_file_file,"r") as infile:
		for line in infile:
			if len(line) < 3: continue
			spls = line.strip().split("\t")
			if spls[0] == "IN": INGROUPS.append(spls[1])
			elif spls[0] == "OUT": OUTGROUPS.append(spls[1])
			else:
				print "Check taxon_code_file file format"
				sys.exit()
	if len(set(INGROUPS) & set(OUTGROUPS)) > 0:
		print "Taxon ID",set(INGROUPS) & set(OUTGROUPS),"in both ingroups and outgroups"
		sys.exit(0)
	print len(INGROUPS),"ingroup taxa and",len(OUTGROUPS),"outgroup taxa read"
	print "Ingroups:",INGROUPS
	print "Outgroups:",OUTGROUPS
	
	for treefile in os.listdir(homoDIR):
		if not treefile.endswith(tree_file_eneding): continue
		with open(homoDIR+treefile,"r") as infile:
			 intree = newick3.parse(infile.readline())
		curroot = intree
		all_names = tree_utils.get_front_names(curroot)
		num_tips = len(all_names)
		num_taxa = len(set(all_names))
		print treefile
		
		#check taxonIDs
		ingroup_names = []
		outgroup_names = []
		for name in all_names:
			if name in INGROUPS:
				ingroup_names.append(name)
			elif name in OUTGROUPS:
				outgroup_names.append(name)
			else:
				print name,"not in ingroups or outgroups"
				sys.exit()
		if len(set(ingroup_names)) < min_ingroup_taxa:
			print "not enough ingroup taxa in tree"
			continue
		
		outID = outDIR + tree_utils.get_clusterID(treefile)
		if len(outgroup_names) > 0: #at least one outgroup present, root and cut inclades
			inclades = tree_utils.extract_rooted_ingroup_clades(curroot,\
				INGROUPS,OUTGROUPS,min_ingroup_taxa)
			inclade_count = 0
			for inclade in inclades:
				inclade_count += 1
				inclade_name = outID+".inclade"+str(inclade_count)
				with open(inclade_name,"w") as outfile:
					outfile.write(newick3.tostring(inclade)+";\n")
				orthologs = tree_utils.get_ortho_from_rooted_inclade(inclade)
				ortho_count = 0
				for ortho in orthologs:
					if len(tree_utils.get_front_labels(ortho)) >= min_ingroup_taxa:
						ortho_count += 1
						with open(inclade_name+".ortho"+str(ortho_count)+".tre","w") as outfile:
							outfile.write(newick3.tostring(ortho)+";\n")

		elif len(all_names) == num_taxa:
			#only output ortho tree when there is no taxon repeats
			with open(outID+".unrooted-ortho.tre","w") as outfile:
				outfile.write(newick3.tostring(curroot)+";\n")
		
		else: #do not attempt to infer direction of gene duplication without outgroup info
			print "duplicated taxa in unrooted tree"
		

if __name__ == "__main__":
	if len(sys.argv) != 6:
		print "python prune_paralogs_RT.py homoTreeDIR tree_file_eneding outDIR min_ingroup_taxa taxon_code_file"
		sys.exit(0)
	
	homoTreeDIR,tree_file_ending,outDIR,min_ingroup_taxa,taxon_code_file=sys.argv[1:]
	RT(homoTreeDIR,tree_file_ending,outDIR,min_ingroup_taxa,taxon_code_file)

