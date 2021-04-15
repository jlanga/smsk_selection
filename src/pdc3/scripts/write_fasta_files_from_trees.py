"""
Read a concatenated fasta file
write individual fasta files for each tree
"""

import sys,os,newick3,phylo3
import tree_utils
from seq import read_fasta_file

def main(fasta,treDIR,tree_file_ending,outDIR):
	if treDIR[-1] != "/": treDIR += "/"
	if outDIR[-1] != "/": outDIR += "/"
	print("Reading fasta file",fasta)
	seqDICT = {} #key is seqID, value is seq
	for s in read_fasta_file(fasta):
		seqDICT[s.name] = s.seq
	print("Writing fasta files")
	filecount = 0
	for i in os.listdir(treDIR):
		if i.endswith(tree_file_ending):
			print(i)
			filecount += 1
			with open(treDIR+i,"r")as infile:
				intree = newick3.parse(infile.readline())
			clusterID = tree_utils.get_clusterID(i)
			if clusterID.endswith("rr"):
				outname = outDIR+clusterID+"_rr.fa"
			else: outname = outDIR+clusterID+"rr.fa"
			with open(outname,"w") as outfile:
				for label in tree_utils.get_front_labels(intree):
					outfile.write(">"+label+"\n"+seqDICT[label]+"\n")
	assert filecount > 0,\
		"No file ends with "+tree_file_ending+" found in "+treDIR

if __name__ =="__main__":
	if len(sys.argv) != 5:
		print("usage: python write_fasta_files_from_trees.py fasta treDIR tree_file_ending outDIR")
		sys.exit()
	
	fasta,treDIR,tree_file_ending,outDIR = sys.argv[1:]
	main(fasta,treDIR,tree_file_ending,outDIR)

