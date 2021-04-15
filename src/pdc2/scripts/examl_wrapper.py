"""
Input: a dir of alignments in phylip format and end with ".phy"
"""

ALIGNMENT_FILE_ENDING = ".phy"

import os,sys

if __name__ == "__main__":
	if len(sys.argv) != 6:
		print "python examl_wrapper.py .phy .model outname number_cores DNA/aa"
		sys.exit(0)
	
	phy = sys.argv[1]
	model = sys.argv[2]
	outname = sys.argv[3]
	num_cores = sys.argv[4]
	seqtype = sys.argv[5]
	
	if seqtype == "aa":
		cmd = "PROT"
	elif seqtype == "DNA":
		m = "DNA" 	
	else:
		print "Input data type: DNA or aa"
		sys.exit()
	
	#build binary input file
	cmd = "parse-examl -s "+phy+" -q "+model+" -m "+m+" -n "+outname
	print cmd
	#os.system(cmd)
	
	#compute a parsimony start tree using raxml
	cmd = "raxmlHPC-PTHREADS-SSE3 -T "+num_cores+" -p 12345 -s "+phy+" -q "+model+" -n "+outname+" -y -m "
	if seqtype == "aa": cmd += "PROTCATWAG"
	else: cmd += "GTRCAT"
	print cmd
	os.system(cmd)
	
	#run examl
	cmd = "mpirun.openmpi -np "+num_cores+" examl -s "+outname+".binary -t RAxML_parsimonyTree."+outname+" -m PSR -n "+outname
	print cmd
	os.system(cmd)

