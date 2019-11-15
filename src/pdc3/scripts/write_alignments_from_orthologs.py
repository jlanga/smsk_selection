import sys,os,newick3,phylo3
import tree_utils
from seq import read_fasta_file

def ortho_to_aln(alndir,tredir,outdir,ortho_tree_file_ending=".tre"):
	"""
	Read final homolog
	write individual alignment files for each ortholog
	Shorten seq id to taxon id
	"""
	if alndir[-1] != "/": alndir += "/"
	if tredir[-1] != "/": tredir += "/"
	if outdir[-1] != "/": outdir += "/"
	filecount = 0
	for i in os.listdir(tredir):
		if i.endswith(ortho_tree_file_ending):
			filecount += 1
			print(i)		
			#read in the alignment into an dictionary
			seqDICT = {} #key is seqID, value is seq
			for s in read_fasta_file(alndir+i.split(".")[0]+".fa.mafft.aln"):
				seqDICT[s.name] = s.seq

			#read in tree tips and write output alignment
			with open(tredir+i,"r")as infile:
				intree = newick3.parse(infile.readline())
			labels = tree_utils.get_front_labels(intree)
			with open(outdir+i.replace(ortho_tree_file_ending,".aln"),"w") as outfile:
				for lab in labels:
					outfile.write(">"+tree_utils.get_name(lab)+"\n"+seqDICT[lab]+"\n")
	assert filecount > 0,\
		"No file ends with "+ortho_tree_file_ending+" was found in "+tredir
		
if __name__ =="__main__":
	if len(sys.argv) != 4:
		print("usage: python write_alignments_from_orthologs.py alnDIR orthoTreeDIR outDIR")
		sys.exit()
	
	alndir,tredir,outdir=sys.argv[1:]
	ortho_to_aln(alndir,tredir,outdir)
