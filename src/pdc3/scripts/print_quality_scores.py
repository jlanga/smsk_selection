"""
Print quality scores to plot in R
"""

import sys,os
import seq
import gzip

SAMPLE_FREQ = 1000

def print_quality_scores(fq):
	infile = gzip.open(fq,"rb") if fq.endswith(".gz") else open(fq,"r")
	outfile = open(fq+".qual","w")
	count = 0
	for read_obj in seq.fastq_generator(infile):
		count += 1 #keep track of number of input read pairs	
		if count % SAMPLE_FREQ == 0:
			outfile.write(" ".join([str(score) for score in read_obj.qualarr])+"\n")
	outfile.close()
	print("quality scores written to "+fq+".qual")
	
if __name__ == "__main__":
	if len(sys.argv) != 2:
		print("Usage: print_quality_scores.py fastq")
		sys.exit(0)
	
	print_quality_scores(sys.argv[1])
	


