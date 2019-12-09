"""
Input:
- prepare a input fasta file that contains sequences used as bait:
	genename.pep.fa
- prepare a directory with candidate sequence pep fasta files named
	taxonid.pep.fa, or
	taxonid.pep.fa.cdhitent
  for each sequence, the sequence ids look like taxonid@seqid

swipe and/or blastp, taking the top 20 hits, and construct the homolog tree
"""

import phylo3,newick3,os,sys
import seq
import ntpath

RUN_BLASTP = True # run blastp search along side of swipe for comparison

def get_filename_from_path(path):
	head, tail = ntpath.split(path)
	return tail or ntpath.basename(head)
	
def swipe(query_fasta,DIR,num_cores,max_num_hits=20,min_bitscore=20.0):
	"""
	given a DIR with peptide fasta files that either end with .pep.fa or 
	.cdhit, swipe on each one
	return a fasta file with hits from each taxa plus the queries added
	"""
	if DIR[-1] != "/": DIR += "/"
	max_num_hits = int(max_num_hits)
	min_bitscore = float(min_bitscore)
	
	# swipe with each taxon
	pepfiles = [i for i in os.listdir(DIR) if (i.endswith(".pep.fa") or i.endswith(".cdhit"))]
	datasets = [i.split(".")[0] for i in pepfiles]
	print len(pepfiles),"input peptide files read"
	assert len(set(datasets)) == len(datasets),\
		"dataset name repeats. remove duplicated sets"
	for i in os.listdir(DIR):
		if i.endswith(".pep.fa") or i.endswith(".cdhit"):
			if not os.path.exists(DIR+i+".psd"):
				os.system("makeblastdb -in "+DIR+i+" -parse_seqids -dbtype prot -out "+DIR+i)
			swipe_outname = DIR+i+"."+get_filename_from_path(query_fasta).split(".")[0]+".swipe"
			if not os.path.exists(swipe_outname):
				cmd = "swipe -d "+DIR+i+" -i "+query_fasta+" -a "+str(num_cores)
				cmd += " -p blastp -o "+swipe_outname+" -m 8 -e 10"
				print cmd
				os.system(cmd)
			assert os.path.exists(swipe_outname), \
				"swipe did not finish correctly"
			"""
			swipe output colums are:
			Query id, Subject id, % identity, alignment length, mismatches, 
			gap openings, q. start, q. end, s. start, s. end, e-value, bit score
			"""
			# summarize the hit seq ids
			if not os.path.exists(swipe_outname+".hits"):
				hit_tuples = [] # alist of tuples (hit, bitscore)
				with open(swipe_outname,"r") as infile:
					for line in infile:
						if len(line) < 3: continue # skip empty lines
						spls = line.strip().split("\t")
						query,hit,bitscore = spls[0],spls[1].replace("lcl|",""),float(spls[-1])
						if query != hit and bitscore >= min_bitscore:
							hit_tuples.append((hit,bitscore))
				
				out = [] # unique hit ids
				for hit, bitscore in sorted(hit_tuples,key=lambda x:x[1],reverse=True):
					if hit not in out:
						out.append(hit)
					if len(out) == max_num_hits:
						break
				if len(out) == 0: print "Warning: No hits found"
				with open(swipe_outname+".hits","w") as outfile:
					for hit in out:
						print hit
						outfile.write(hit+"\n")
	
	# write output fasta
	outname = query_fasta.replace(".pep.fa","_swipe.fa")
	print "Writing output fasta",outname
	outfile = open(outname,"w")
	query_seqids = [] # avoid seq id repeats
	with open(query_fasta,"r") as infile:
		for line in infile:
			outfile.write(line) # copy over query seqs
			if line[0] == ">":
				query_seqids.append(line.strip()[1:])
	for i in os.listdir(DIR):
		if i.endswith(".pep.fa") or i.endswith(".cdhit"):
			seqDICT = {} # key is seq name, value is seq
			for s in seq.read_fasta_file(DIR+i):
				seqDICT[s.name] = s.seq
			with open(DIR+i+"."+get_filename_from_path(query_fasta).split(".")[0]+".swipe.hits","r") as infile:
				for line in infile:
					line = line.strip()
					if len(line) > 0 and line not in query_seqids:
						outfile.write(">"+line+"\n"+seqDICT[line]+"\n")
	outfile.close()

def blastp(query_fasta,DIR,num_cores,max_num_hits=20,min_bitscore=20.0):
	"""
	same as swipe but using blastp
	"""
	if DIR[-1] != "/": DIR += "/"
	max_num_hits = int(max_num_hits)
	min_bitscore = float(min_bitscore)
	
	# blastp with each taxon
	pepfiles = [i for i in os.listdir(DIR) if (i.endswith(".pep.fa") or i.endswith(".cdhit"))]
	datasets = [i.split(".")[0] for i in pepfiles]
	print len(pepfiles),"input peptide files read"
	assert len(set(datasets)) == len(datasets),\
		"dataset name repeats. remove duplicated sets"
	for i in os.listdir(DIR):
		if i.endswith(".pep.fa") or i.endswith(".cdhit"):
			if not os.path.exists(DIR+i+".psd"):
				os.system("makeblastdb -in "+DIR+i+" -parse_seqids -dbtype prot -out "+DIR+i)
			blastp_outname = DIR+i+"."+get_filename_from_path(query_fasta).split(".")[0]+".blastp"
			if not os.path.exists(blastp_outname):
				cmd = "blastp -db "+DIR+i
				cmd += " -query "+query_fasta
				cmd += " -num_threads "+str(num_cores)
				cmd += " -out "+blastp_outname
				cmd += " -evalue 10"
				cmd += " -max_target_seqs "+str(max_num_hits)
				cmd += " -outfmt '6 qseqid qlen sseqid slen frames pident nident length mismatch gapopen qstart qend sstart send evalue bitscore'"
				print cmd
				os.system(cmd)
			assert os.path.exists(blastp_outname), \
				"blastp did not finish correctly"
			"""
			blastp output colums are:
			0-qseqid 2-sseqid 15-bitscore'"
			"""
			# summarize the hit seq ids
			if not os.path.exists(blastp_outname+".hits"):
				hit_tuples = [] # alist of tuples (hit, bitscore)
				with open(blastp_outname,"r") as infile:
					for line in infile:
						if len(line) < 3: continue # skip empty lines
						spls = line.strip().split("\t")
						query,hit,bitscore = spls[0],spls[2],float(spls[-1])
						if query != hit and bitscore >= min_bitscore:
							hit_tuples.append((hit,bitscore))
				
				out = [] # unique hit ids
				for hit, bitscore in sorted(hit_tuples,key=lambda x:x[1],reverse=True):
					if hit not in out:
						out.append(hit)
					if len(out) == max_num_hits:
						break
				if len(out) == 0: print "Warning: No hits found"
				with open(blastp_outname+".hits","w") as outfile:
					for hit in out:
						print hit
						outfile.write(hit+"\n")
	
	# write output fasta
	outname = query_fasta.replace(".pep.fa","_blastp.fa")
	print "Writing output fasta",outname
	outfile = open(outname,"w")
	query_seqids = [] # avoid seq id repeats
	with open(query_fasta,"r") as infile:
		for line in infile:
			outfile.write(line) # copy over query seqs
			if line[0] == ">":
				query_seqids.append(line.strip()[1:])
	for i in os.listdir(DIR):
		if i.endswith(".pep.fa") or i.endswith(".cdhit"):
			seqDICT = {} # key is seq name, value is seq
			for s in seq.read_fasta_file(DIR+i):
				seqDICT[s.name] = s.seq
			with open(DIR+i+"."+get_filename_from_path(query_fasta).split(".")[0]+".blastp.hits","r") as infile:
				for line in infile:
					line = line.strip()
					if len(line) > 0 and line not in query_seqids:
						outfile.write(">"+line+"\n"+seqDICT[line]+"\n")
	outfile.close()

if __name__ == "__main__":
	if len(sys.argv) != 5:
		print "python bait_homologs.py query_pep_fa databaseDIR num_to_bait num_cores"
		sys.exit(0)
	
	query_fasta,DIR,num_to_bait,num_cores = sys.argv[1:]
	swipe(query_fasta,DIR,num_cores,max_num_hits=num_to_bait)
	if RUN_BLASTP: blastp(query_fasta,DIR,num_cores,max_num_hits=num_to_bait)
