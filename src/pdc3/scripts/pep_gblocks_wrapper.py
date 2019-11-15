"""
assuming that all alignment files end with ".aln"
"""

import os,sys

FILE_ENDING = ".aln"

#seqtype = "--nuc"
seqtype = "--amino"

MIN_CHR = 40

if __name__ == "__main__":
	if len(sys.argv) != 3:
		print("usage: python pep_gblocks_wrapper.py inDIR outDIR")
		sys.exit()
	
	inDIR = sys.argv[1]+"/"
	outDIR = sys.argv[2]+"/"
	for i in os.listdir(inDIR):
		if i[-len(FILE_ENDING):] != FILE_ENDING: continue
		
		#count how many seqs in the alignment
		seq_count = 0
		with open(inDIR+i,"r") as infile:
			for line in infile:
				if line[0] == ">": seq_count += 1
		min_col = seq_count/2 + 1 #bigger than half
		
		#b1 min number of chr for a conserved position
		#b2 min number of chr for a flank position
		#b3 max number of contiguous nonconserved positions
		#b4 min length of a block
		#b5 allowed gap positions h=with half, a=all
		#p results and parameter files n=no
		com = "Gblocks "+inDIR+i+" -t=p -b1="+str(min_col)+" -b2="+str(min_col)
		com += " -b3=8 -b4=10 -b5=h -p=n" #-t=d for DNA
		print(com)
		os.system(com)
			
		#remove empty columns produced by Gblocks
		cmd = "phyutility -aa -clean 0.00001 -in "+inDIR+i
		cmd += "-gb -out "+outDIR+i+"-cln"
		print(cmd)
		os.system(cmd)
		
		#remove intermediate files
		os.system("rm "+inDIR+i+"-gb")
