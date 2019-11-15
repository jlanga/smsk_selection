"""
make sure that "raxml" is in the path

Input: a dir of cleaned alignments in fasta format and end with ".aln-cln"
Output: tree names are clusterID.raxml.tre
"""

import os,sys
import subprocess
from seq import read_fasta_file

def raxml(DIR,cleaned,num_cores,seqtype):
	assert cleaned.endswith(".aln-cln"),\
		"raxml infile "+cleaned+" not ends with .aln-cln"
	assert seqtype == "aa" or seqtype == "dna","Input data type: dna or aa"
	if len(read_fasta_file(DIR+cleaned)) < 4:
		sys.stderr.write("less than 4 sequences in "+DIR+cleaned)
		return
	clusterID = cleaned.split(".")[0]
	tree = DIR+clusterID+".raxml.tre"
	raw_tree = "RAxML_bestTree."+cleaned
	model = "PROTCATWAG" if seqtype == "aa" else "GTRCAT"
	if not os.path.exists(tree) and not os.path.exists(raw_tree):
		# raxml crashes if input file starts with . 
		infasta = cleaned if DIR == "./" else DIR+cleaned
		cmd = ["raxml","-T",str(num_cores),"-p","12345","-s",\
			   infasta,"-n",cleaned,"-m",model]
		print(" ".join(cmd))
		p = subprocess.Popen(cmd,stdout=subprocess.PIPE)
		out = p.communicate()
		assert p.returncode == 0,"Error raxml"+out[0]
	try:
		os.rename(raw_tree,tree)
		os.remove("RAxML_info."+cleaned)
		os.remove("RAxML_log."+cleaned)
		os.remove("RAxML_parsimonyTree."+cleaned)
		os.remove("RAxML_result."+cleaned)
		os.remove(DIR+cleaned+".reduced")
	except: pass # no need to worry about extra intermediate files
	return tree

def main(DIR,num_cores,seqtype):
	if DIR[-1] != "/": DIR += "/"
	filecount = 0
	for i in os.listdir(DIR):
		if i.endswith(".aln-cln"):
			filecount += 1
			raxml(DIR,i,num_cores,seqtype)
	assert filecount > 0, "No file end with .aln-cln found in "+DIR
	
	
if __name__ == "__main__":
	if len(sys.argv) != 4:
		print("python raxml_wrapper.py DIR number_cores dna/aa")
		print("make sure that the executable is named 'raxml' and is in the path")
		sys.exit(0)
	
	DIR, num_cores, seqtype  = sys.argv[1:]
	main(DIR,num_cores,seqtype)
	
