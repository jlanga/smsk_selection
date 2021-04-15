"""
Filter the aligned and cleaned ortholog matrices by number of taxa and characters
Write out name of matrices that passed the filter
Also write out supermatrix stats
"""
from seq import read_fasta_file
from tree_utils import get_name
import sys,os

def concatenate(clnDIR,numofsitesFilter,numoftaxaFilter,outname):
	"""filter cleaned alignments and concatenate"""
	if clnDIR[-1] != "/": clnDIR += "/"
	sites_filter = int(numofsitesFilter)
	taxa_filter = int(numoftaxaFilter)
	
	print "Filtering ortholog matrixes"
	selected = [] # list of alignment file names that pass the filters
	for i in os.listdir(clnDIR):
		if i.endswith(".aln-cln"):
			seqlist = read_fasta_file(clnDIR+i)
			num_seq = len(seqlist)
			num_col = len(seqlist[0].seq)
			if num_seq >= taxa_filter and num_col >= sites_filter:
				selected.append(i)
	print len(selected),"matrices passed the filter"

	print "Getting matrix occupancy stats"
	taxon_occupancy = {}
	#key is taxon name, value is [times present in a matrix,total length for this taxon]
	total_aligned_len = 0 #record how long the final concatenated matrix is
	
	cmd = "pxcat"+" -o "+outname+".fa"+" -p "+outname+".model"+" -s "
	for i in selected:
		cmd += clnDIR+i+" "
		seqlist = read_fasta_file(clnDIR+i)
		total_aligned_len += len(seqlist[0].seq)
		for s in seqlist:
			taxonid = get_name(s.name)
			if taxonid not in taxon_occupancy:
				taxon_occupancy[taxonid] = [0,0]
			taxon_occupancy[taxonid][0] += 1
			taxon_occupancy[taxonid][1] += len((s.seq.replace("-","")).replace("?",""))
	cmd += "\n"
	
	total_ortho = len(selected)
	with open(outname+"_taxon_occupancy_stats","w") as outfile:
		outfile.write("taxon\t#orthologs\t#total_charactors\tperc_orthologs\tperc_charactors\n")
		sum_char = 0
		for taxon in taxon_occupancy:
			times,chars = taxon_occupancy[taxon][0],taxon_occupancy[taxon][1]
			sum_char += chars
			out = taxon+"\t"+str(times)+"\t"+str(chars)+"\t"
			out += str(times/float(total_ortho))+"\t"+str(chars/float(total_aligned_len))+"\n"
			outfile.write(out)
		total_taxa = len(taxon_occupancy)
		out = "\nSupermatrix dimension "+str(total_taxa)+" taxa, "
		out += str(total_ortho)+" loci and "+str(total_aligned_len)+" aligned columns\n"
		out += "Overall matrix occupancy "+str(sum_char/float(total_taxa*total_aligned_len))+"\n"
		outfile.write(out)

	print "Supermatrix taxon occupancy stats written to",outname+"_taxon_occupancy_stats"
	print "Waiting for concatenation to finish. This may take several minutes..."
	with open(outname+".temp.sh","w") as f: f.write(cmd)
	os.system("bash "+outname+".temp.sh")	
	
	#writes phy file
	cmd_pxs2phy = ["pxs2phy","-o",outname+".phy","-s",outname+".fa"]
	print (" ".join(cmd_pxs2phy))
	os.system(" ".join(cmd_pxs2phy))
	
	#writes nex file
	cmd_pxs2nex = ["pxs2nex","-o",outname+".nex","-s",outname+".fa"]
	print (" ".join(cmd_pxs2nex))
	os.system(" ".join(cmd_pxs2nex))

	
	assert os.path.exists(outname+".phy") and os.path.exists(outname+".nex") and os.path.exists(outname+".fa"),  "error concatenate"
	os.system("rm "+outname+".temp.sh")
	print "outfiles written",outname+".phy",outname+".model"

if __name__ == "__main__":
	if len(sys.argv) != 5:
		print "usage: python concatenate_matrices_phix.py aln-clnDIR numofsitesFilter numoftaxaFilter outname"
		sys.exit()
		
	clnDIR,numofsitesFilter,numoftaxaFilter,outname = sys.argv[1:]
	concatenate(clnDIR,numofsitesFilter,numoftaxaFilter,outname)
