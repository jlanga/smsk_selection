import sys,os

"""
input: a DIR with fasta files of peptide sequences
output: all-by-all blastp results in customized tabular format
"""

if __name__ == "__main__":
	if len(sys.argv) != 4:
		print "python all-by-all_blastp.py DIR file_ending num_threads"
		sys.exit(0)
	
	DIR = sys.argv[1]
	if DIR[-1] != "/": DIR = DIR+"/"
	FILE_ENDING = sys.argv[2]
	NUM_CORES = sys.argv[3]
	
	#concatenate all the fasta files
	cmd = "cat "+DIR+"*"+FILE_ENDING+" >"+DIR+"all.fa"
	print cmd
	os.system(cmd)

	#create blast database
	print "Creating blast database"
	cmd = "/home/yangya/dmorales/apps/ncbi-blast-2.7.1+/bin/makeblastdb -in "+DIR+"all.fa -parse_seqids -dbtype prot -out "+DIR+"all.fa"
	print cmd
	os.system(cmd)

	#blast against the database one by one
	#setting the max_target_seqs to 1000 to avoid braking up large gene families in MCL
	#setting the evlalue cutoff to 10 to get as many hits as possible
	for i in os.listdir(DIR):
		if i.endswith(FILE_ENDING) and i[:3] != "all":
			print "conducting blast on "+i
			cmd = "/home/yangya/dmorales/apps/ncbi-blast-2.7.1+/bin/blastp -db "+DIR+"all.fa -query "+DIR+i
			cmd += " -evalue 10 -num_threads "+NUM_CORES
			cmd += " -max_target_seqs 1000 -out "+DIR+i+".rawblastp "
			cmd += "-outfmt '6 qseqid qlen sseqid slen frames pident nident length mismatch gapopen qstart qend sstart send evalue bitscore'"
			print cmd
			os.system(cmd)
	print "blast is done"

