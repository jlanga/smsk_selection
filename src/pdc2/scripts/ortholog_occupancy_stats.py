import sys,os,newick3,phylo3

"""
Get ortholog gene occupancy stats
"""

file_ending = ".tre"

def get_name(label):
	if "@" in label:
		return label.split("@")[0]
	else: return ""
	
def get_front_labels(node):
	leaves = node.leaves()
	return [i.label for i in leaves]
	
def get_front_names(node):
	labels = get_front_labels(node)
	return [get_name(i) for i in labels]
	
if __name__ == "__main__":
	if len(sys.argv) != 2:
		print "usage: python ortholog_occupancy_stats.py ortho_treDIR"
		sys.exit(0)
	DIR = sys.argv[1]+"/"
	outfile = open("ortho_stats","w")
	DICT = {} #key is taxon name, value is how many orthologs it is in
	total_ortho = 0
	for i in os.listdir(DIR):
		if i[-len(file_ending):] == file_ending and "ortho" in i:
			print i
			total_ortho += 1
			with open(DIR+i,"r") as infile:
				intree = newick3.parse(infile.readline())
			names = get_front_names(intree)
			for taxon in names:
				if taxon not in DICT:
					DICT[taxon] = 0
				DICT[taxon] += 1
			outfile.write(str(len(names))+"\n")
	outfile.close()
	print "number of taxa in each ortholog written to ortho_stats"
	
	with open("taxon_stats","w") as outfile:
		outfile.write("taxonID\tnum_ortho\t%ortho_out_of_total_"+str(total_ortho)+"\n")
		for taxon in DICT:
			outfile.write(taxon+"\t"+str(DICT[taxon])+"\t"+str(float(DICT[taxon])/total_ortho)+"\n")
	print "number of ortholog for each taxon written to taxon_stats"
			
