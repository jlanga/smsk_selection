"""
input: filtered blast hits filtered by hit fraction
output connected components
"""

import os,sys
import networkx as nx

def get_taxonID(seqID):
	return seqID.split("@")[0]
	
if __name__ == "__main__":
	if len(sys.argv) != 3:
		print("usage: python find_connected_components.py edges minimal_taxa")
		sys.exit()
	
	edges = sys.argv[1]
	MIN_TAXA = int(sys.argv[2])
	DIR = sys.argv[3]
	
	print("Reading edges")
	infile = open(edges,"r")
	G = nx.Graph() #initiate the graph
	for line in infile:
		if len(line) > 3:
			spls = line.strip().split("\t")
			G.add_edge(spls[0],spls[1])
	infile.close()
	
	print("Writing connected components")
	ccs = nx.connected_components(G) #get the connected components
	outfile1 = open(edges+".ccs","w")
	
	for cc in ccs: #cc is a list of seqids in
		SG = G.subgraph(cc)
		#check how many taxa are in the connected component
		taxa = [] #store a list of unique taxon codes
		for seqid in cc:
			taxonID = get_taxonID(seqid)
			if taxonID not in taxa:
				taxa.append(taxonID)
		if len(taxa) >= MIN_TAXA:
			for seqid in cc:
				outfile1.write(seqid+"\t")
			outfile1.write("\n")
	outfile1.close()


