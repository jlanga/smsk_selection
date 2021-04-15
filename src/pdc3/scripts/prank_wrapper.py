import os,sys

if __name__ == "__main__":
	if len(sys.argv) != 5:
		print("usage: python prank_wrapper.py inDIR outDIR infile_ending dna/codon/aa")
		sys.exit()
	
	inDIR = sys.argv[1]+"/"
	outDIR = sys.argv[2]+"/"
	file_end = sys.argv[3]
	
	if sys.argv[4] == "aa":
		seqtype = "-protein"
	elif sys.argv[4] == "codon":
		seqtype = "-codon"
	elif sys.argv[4] == "dna":
		seqtype = "-DNA"
	else:
		print("Input data type: DNA or aa")
		sys.exit()
	
	done = os.listdir(outDIR)
	for i in os.listdir(inDIR):
		if i[-len(file_end):] != file_end: continue
		if i+".aln" in done: continue
		
		#replace U from the sequences to avoid crashing prank
		infile = open(inDIR+i,"r")
		outfile = open(inDIR+i+".temp","w")
		for line in infile:
			if line[0] != ">":
				line = line.replace("U","X")
				line = line.replace("u","x")
			outfile.write(line)
		outfile.close()
		infile.close()
		
		out = outDIR+i+".aln "
		cmd = "prank -d="+inDIR+i+".temp -o="+out+" "+seqtype
		print(cmd)
		os.system(cmd)
		os.system("rm "+inDIR+i+".temp")
		os.system("mv "+outDIR+i+".aln.best.fas "+out)
