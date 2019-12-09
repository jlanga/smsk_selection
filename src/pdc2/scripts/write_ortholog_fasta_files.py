"""
Read initial concatenated fasta before clustering

Since no taxon repeats
Shorten seq id to taxon id
"""

import sys,os,newick3,phylo3
from Bio import SeqIO

ORTHO_TREE_FILE_ENDING = ".tre"

def get_name(label):
	return label.split("@")[0]

def get_front_labels(node):
	leaves = node.leaves()
	return [i.label for i in leaves]

if __name__ =="__main__":
	if len(sys.argv) != 5:
		print "usage: python write_ortholog_fasta_files.py fasta treeDIR outDIR MIN_TAXA"
		sys.exit()
	
	fasta = sys.argv[1]
	treDIR = sys.argv[2]+"/"
	outDIR = sys.argv[3]+"/"
	MIN_TAXA = int(sys.argv[4])

	print "Reading the original fasta file"
	handle = open(fasta,"rU")
	#hash table of taxonID -> seqID -> seq
	seqDICT = {} #key is taxonID, value is seqID
	for seq_record in SeqIO.parse(handle,"fasta"):
		seqID,seq = str(seq_record.id),str(seq_record.seq)
		taxonID = get_name(seqID)
		if taxonID not in seqDICT:
			seqDICT[taxonID] = {} #key is taxonID, value is seq
			print "Adding sequence from",taxonID
		seqDICT[taxonID][seqID] = seq
	handle.close()
	
	print "Writing fasta files"
	for i in os.listdir(treDIR):
		if i[-len(ORTHO_TREE_FILE_ENDING):] == ORTHO_TREE_FILE_ENDING:
			#read in tree tips and write output alignment
			with open(treDIR+i,"r")as infile:
				intree = newick3.parse(infile.readline())
			labels = get_front_labels(intree)
			if len(labels) >= MIN_TAXA:
				with open(outDIR+i.replace(ORTHO_TREE_FILE_ENDING,".fa"),"w") as outfile:
					for lab in labels:
						name = get_name(lab)
						outfile.write(">"+name+"\n"+seqDICT[name][lab]+"\n")

