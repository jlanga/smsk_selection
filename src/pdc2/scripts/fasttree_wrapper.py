import os,sys
import subprocess

def fasttree(DIR,cleaned,seqtype):
	"""read in a cleaned alignment ends with '.aln-cln
	extimate a tree using fasttree'"""
	if DIR[-1] != "/": DIR += "/"
	assert cleaned.endswith(".aln-cln"),\
		"fasttree infile "+cleaned+" not ends with .aln-cln"
	assert seqtype == "aa" or seqtype == "dna","Input data type: dna or aa"
	tree = cleaned.split(".")[0]+".fasttree.tre"
	alg = ["-wag"] if seqtype == "aa" else ["-nt","-gtr"]
	if os.path.exists(DIR+tree):
		return DIR+tree
	cmd = ["fasttree"]+alg+["-quiet",DIR+cleaned]
	out = open(DIR+tree, 'w')
	p = subprocess.Popen(cmd,stdout=out)
	out.close()
	p.communicate()
	assert p.returncode == 0,"Error fasttree"
	return DIR+tree

if __name__ == "__main__":
	if len(sys.argv) != 3:
		print "python fasttree_wrapper.py DIR dna/aa"
		print "make sure the executable for fasttree is in the path and is called 'fasttree'"
		sys.exit(0)

	DIR, seqtype  = sys.argv[1:]
	if DIR[-1] != "/": DIR += "/"
	filecount = 0
	for i in os.listdir(DIR):
		if i.endswith(".aln-cln"):
			filecount += 1
			fasttree(DIR,i,seqtype)
	assert filecount > 0, "No file end with .aln-cln found in "+DIR
				
