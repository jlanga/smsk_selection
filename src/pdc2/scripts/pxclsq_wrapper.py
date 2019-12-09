"""
pxclsq (from phyx) wrapper. adapted from Ya's phyutility_wrapper.py
make sure that pxclsq is in the path
"""

import sys,os
from seq import read_fasta_file

def pxclsq(DIR,alignment,min_col_occup,seqtype,min_chr=10):
	"""
	remove columns with occupancy lower than MIN_COLUMN_OCCUPANCY
	remove seqs shorter than MIN_CHR after filter columns
	"""
	if DIR[-1] != "/": DIR += "/"
	cleaned = alignment+"-cln"
	if os.path.exists(DIR+cleaned): return cleaned
	assert alignment.endswith(".aln"),\
		"pxclsq infile "+alignment+" not ends with .aln"
	assert os.stat(DIR+alignment).st_size > 0, DIR+alignment+"empty"
	assert seqtype == "aa" or seqtype == "dna","Input data type: dna or aa"

	if seqtype == "aa":
		cmd = ["pxclsq","-a","-p",str(min_col_occup),"-s",\
			   DIR+alignment,"-o",DIR+alignment+"-pht"]
	else:
		cmd = ["pxclsq","-p",str(min_col_occup),"-s",\
			   DIR+alignment,"-o",DIR+alignment+"-pht"]
	print " ".join(cmd)
	os.system(" ".join(cmd))
	assert os.path.exists(DIR+alignment+"-pht"),"Error pxclsq"
	
	#remove empty and very short seqs
	outfile = open(DIR+cleaned,"w")
	for s in read_fasta_file(DIR+alignment+"-pht"):
		if len(s.seq.replace("-","")) >= min_chr:
			outfile.write(s.get_fasta())
	outfile.close()
	os.remove(DIR+alignment+"-pht")
	return cleaned
	
def main(DIR,min_col_occup,seqtype):
	if DIR[-1] != "/": DIR += "/"
	filecount = 0
	for i in os.listdir(DIR):
		if i.endswith(".aln"):
			filecount += 1
			pxclsq(DIR,i,min_col_occup,seqtype)
	assert filecount > 0, "No file end with .aln found in "+DIR
				

if __name__ == "__main__":
	if len(sys.argv) != 4:
		print "python pxclsq_wrapper.py DIR min_column_occupancy dna/aa"
		sys.exit(0)

	DIR,min_col_occup,seqtype = sys.argv[1:]
	main(DIR,min_col_occup,seqtype)
