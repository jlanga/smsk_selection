"""
make sure that raxml,  mafft,  run_pasta.py,  pxclsq and fastree are in the path
"""

import os, sys
from mafft_wrapper import mafft
from pasta_wrapper import pasta
from pxclsq_wrapper import pxclsq
from fasttree_wrapper import fasttree
from raxml_wrapper import raxml
from raxml_bs_wrapper import raxml_bs
from seq import read_fasta_file

NUM_SEQ_CUTOFF = 1000 # Use different alignment and tree inference tools
# below vs. above this cutoff

def get_fasta_size(fasta):
	"""
	given a fasta file
	output the number of seqs and the length of the longest seq
	"""
	longest = 0
	seqlist = read_fasta_file(fasta)
	for s in seqlist:
		longest = max(longest, len(s.seq.replace("-", "")))
	return len(seqlist), longest
	
def fasta_to_tree(DIR, fasta, num_cores, seqtype, num_seq_cutoff=NUM_SEQ_CUTOFF):
	"""
	given a fasta file
	align,  trim alignment and build a tree
	choose appropriate tools depending on size of the fasta file
	"""
	if DIR[-1] != "/":
		DIR += "/"
	seqcount,  maxlen = get_fasta_size(DIR+fasta)
	if seqcount < 4:
		sys.stderr.write("Less than four sequences in "+DIR+fasta + "\n")
		return
	print(fasta, seqcount, "sequences")
	if seqcount >= NUM_SEQ_CUTOFF: # large cluster
		print("running pasta")
		alignment = pasta(DIR, fasta, num_cores, seqtype)
		cleaned = pxclsq(DIR, alignment, 0.1, seqtype)
		if len(read_fasta_file(DIR+cleaned)) >= 4:
			tree = fasttree(DIR, cleaned, seqtype)
		else:
			sys.stderr.write("Less than 4 taxa in ", cleaned, "\n")
			return
	else: # small cluster
		alignment = mafft(DIR, fasta, num_cores, seqtype)
		cleaned = pxclsq(DIR, alignment, 0.1, seqtype)
		if len(read_fasta_file(DIR+cleaned)) >= 4:
			tree = raxml(DIR, cleaned, num_cores, seqtype)
		else:
			sys.stderr.write("Less than 4 taxa in", cleaned, "\n")
			return

def fasta_to_bs_tree(DIR, fasta, num_cores, seqtype):
	"""
	given a fasta file for the final homolog
	align,  trim alignment and build a tree with bootstrap support
	"""
	if DIR[-1] != "/": DIR += "/"
	seqcount,  maxlen = get_fasta_size(DIR+fasta)
	if seqcount < 4:
		sys.stderr.write("Less than four sequences in "+DIR+fasta)
		return
	print(fasta, seqcount, "sequences")
	alignment = mafft(DIR, fasta, num_cores, seqtype)
	cleaned = pxclsq(DIR, alignment, 0.2, seqtype)
	if len(read_fasta_file(DIR+cleaned)) >= 4:
		tree = raxml_bs(DIR, cleaned, num_cores, seqtype)
	else: print("Less than 4 taxa in", cleaned)

def main(DIR, num_cores, seqtype, bs, test=False):
	"""if test,  only process clusterID that ends with 0"""
	assert seqtype == "aa" or seqtype == "dna", \
		"Input data type: dna or aa"
	assert bs == "y" or bs=="n", \
		"bootstrap? (y/n)"
	if DIR[-1] != "/": DIR += "/"
	
	#check for really long sequences or large alignments.
	#These crashes the alignment program
	for i in os.listdir(DIR):
		if i.endswith(".fa"):
			seqcount, maxlen = get_fasta_size(DIR+i)
			if (maxlen>=10000 and seqtype=="aa") or (maxlen>=30000 and seqtype=="dna"):
				print(i, "has", seqcount, "sequences")
				print("longest sequence has", maxlen, "characters")
				print("Warning: sequence too long. May crash alignment process")
				#sys.exit()

	filecount = 0
	for i in os.listdir(DIR):
		if not i.endswith(".fa"): continue
		if test and (i.split(".")[0])[-1] != "0": continue
		filecount += 1
		if bs == "n":
			fasta_to_tree(DIR=DIR, fasta=i, num_cores=num_cores, seqtype=seqtype)
		else: fasta_to_bs_tree(DIR=DIR, fasta=i, num_cores=num_cores, seqtype=seqtype)
	assert filecount > 0,  "No file end with .fa found in "+DIR

if __name__ == "__main__":
	if len(sys.argv) != 5:
		print("python fasta_to_tree.py DIR number_cores dna/aa bootstrap(y/n)")
		sys.exit(0)
	
	main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])


