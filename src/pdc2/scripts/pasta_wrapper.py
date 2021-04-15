import sys,os
import subprocess

def pasta(DIR,fasta,thread,seqtype):
	if DIR[-1] != "/": DIR += "/"
	if DIR == "./": DIR = ""
	aln = DIR+fasta+".pasta.aln"
	temp_aln = DIR+"pastajob.marker001."+fasta+".temp.aln"
	fa = DIR+fasta
	if os.path.exists(aln): return fasta+".pasta.aln"
	if not os.path.exists(temp_aln):
		assert seqtype == "aa" or seqtype == "dna","Input data type: dna or aa"
		seq = "Protein" if seqtype == "aa" else "DNA" 	
		
		# pasta does not recognize "*" or "U"
		seqcount = 0 #record how many sequences in the fasta file
		infile = open(fa,"r")
		outfile = open(fa+".temp","w")
		for line in infile:
			if len(line.strip()) == 0: continue # skip empty lines
			if line[0] == ">":
				seqcount += 1
				outfile.write(line)
			else:
				if seqtype == "aa":
					#remove U which is usually not in aa alphabet
					line = (line.replace("U","X")).replace("u","x")
					#remove stop codon and seq after it
					if "*" in line:
						line = line[:line.find("*")]
					if line[-1] != "\n": line += "\n"
				outfile.write(line)
		infile.close()
		outfile.close()
				
		cmd = ["run_pasta.py","--input="+fa+".temp","--datatype="+seq]
		print " ".join(cmd)
		p = subprocess.Popen(cmd)
		err,out = p.communicate()
		print err,out
		assert p.returncode == 0 and os.path.exists(temp_aln),"Error pasta"
	
	#remove intermediate files
	os.remove(fa+".temp")
	os.rename(temp_aln,aln)
	os.system("rm "+DIR+"pastajob*")
	return fasta+".pasta.aln"
	
if __name__ =="__main__":
	if len(sys.argv) != 5:
		print "usage: python pasta_wrapper.py DIR infile_ending num_cores dna/aa"
		sys.exit()
	
	DIR,infile_ending,thread,seqtype = sys.argv[1:]
	if DIR[-1] != "/": DIR += "/"
	filecount = 0
	for i in os.listdir(DIR):
		if i.endswith(infile_ending):
			filecount += 1
			pasta(DIR=DIR,fasta=i,thread=thread,seqtype=seqtype)
	assert filecount > 0, "No file end with "+file_end+"found in "+DIR
				
