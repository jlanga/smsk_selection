"""
Takes a directory of fasta files

If there are >= 1000 sequences in the direction, use --auto
For fasta files with less than 1000 sequences, use the slower but much 
more accurate algorithm

Uncomment the com += "--anysymbol " line if there are "U" or any other unusual
charactors in the sequences
"""

import os,sys
import subprocess
from seq import read_fasta_file

def mafft(DIR,fasta,thread,seqtype):
	if DIR[-1] != "/": DIR += "/"
	alignment = fasta+".mafft.aln"
	if os.path.exists(DIR+alignment) and os.stat(DIR+alignment).st_size>0:
		return alignment
	assert seqtype == "aa" or seqtype == "dna","Input data type: dna or aa"	
	seqlist = read_fasta_file(DIR+fasta)
	seqcount = len(seqlist)
	maxlen = 0
	for s in seqlist:
		maxlen = max(maxlen,len(s.seq))
	assert seqcount >= 2, "less than 4 sequences in "+DIR+fasta
	
	if seqtype == "dna":
		infasta = DIR+fasta
		seq = "--nuc" 
	else:
		infasta = DIR+fasta+".temp"
		seq = "--amino"
		with open(infasta,"w") as outfile:
			for s in seqlist:
				#remove U which is usually not in aa alphabet
				s.seq = s.seq.replace("U","X")
				s.seq = s.seq.replace("u","x")
				#remove stop codon and seq after it
				if "*" in s.seq:
					s.seq = s.seq[:s.seq.find("*")]
				outfile.write(s.get_fasta())

	if seqcount >= 1000 or maxlen >= 10000:
		alg = ["--auto"] #so that the run actually finishes!
	else: alg = ["--genafpair","--maxiterate","1000"]
	
	cmd = ["mafft"]+alg+[seq,"--thread",str(thread)]
	#com += ["--anysymbol"] # when there are "U"s in aa sequences
	cmd += [infasta]
	print " ".join(cmd)
	out = open(DIR+alignment, 'w')
	p = subprocess.Popen(cmd,stderr=subprocess.PIPE,stdout=out)
	out.close()
	p.communicate()
	assert p.returncode == 0,"Error mafft"
	if seqtype == "aa": os.remove(DIR+fasta+".temp")
	return alignment

def main(DIR,infile_ending,thread,seqtype):
	if DIR[-1] != "/": DIR += "/"
	filecount = 0
	for i in os.listdir(DIR):
		if i.endswith(infile_ending):
			filecount += 1
			mafft(DIR=DIR,fasta=i,thread=thread,seqtype=seqtype)
	assert filecount > 0, "No file end with "+file_end+"found in "+DIR
	
if __name__ == "__main__":
	if len(sys.argv) != 5:
		print "usage: python mafft_wrapper.py DIR infile_ending thread dna/aa"
		sys.exit()
	
	DIR,infile_ending,thread,seqtype = sys.argv[1:]
	main(DIR,infile_ending,thread,seqtype)

				
