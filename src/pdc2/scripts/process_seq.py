"""
Modify these following paths according to the location that each file are located:
SRC_HOME
APPS_HOME
BLASTP_DB_PATH
"""

import os,sys
import seq
import gzip
import shutil

APPS_HOME = "/home/yangya/dmorales/apps/" # where trinity and trnasdecoder dirs are located
BLASTP_DB_PATH = APPS_HOME+"/home/yangya/dmorales/data/transdecoder_blastp_db/db" # the custom blast database
TRINITY_CMD = "Trinity"
TRANSDECODER_PFEFIX = APPS_HOME+"TransDecoder-v5.0.2/TransDecoder."
MIN_FASTA = 1000 # miminal number of seqs expected in assembly and translation files

# Pfam settings. Currently disabled
#PFAM_PATH = APPS_HOME+"/data/pfam/Pfam-AB.hmm.bin"
#Setups for bigsmith locally
#HMMSCAN_PATH = ""
# Setups for smithbigmem
#HMMSCAN_PATH = "~/apps/hmmer-3.1b1-linux-intel-x86_64/binaries/hmmscan"


def run(cmd,logfile):
	"""print, run and log calls"""
	print cmd
	os.system(cmd)
	with open(logfile,"a") as outfile: outfile.write(cmd+"\n")

def fasta_ok(fasta,min_count=MIN_FASTA):
	"""count number of non-empty fasta sequences"""
	if not os.path.exists(fasta):
		return False
	fasta_count = 0
	for i in seq.read_fasta_file(fasta):
		if len(i.seq) > 0: fasta_count += 1
	print fasta,"contains",fasta_count,"non-empty sequences"
	if fasta_count >= min_count: return True
	else: return False


def blastpout_ok(blastpout,min_count=MIN_FASTA):
	"""count number of unique query ids"""
	if not os.path.exists(blastpout): return False
	with open(blastpout) as infile:
		count = len(set([line.split("\t")[0] for line in infile]))
	print blastpout,"contains",count,"unique query ids"
	if count >= min_count:
		return True # most fasta should have a hit
	else: return False
	
			
def get_input_string_pe(inDIRs,read1_identifier,read2_identifier,infile_end):
	"""Collect all the pair end fastq files
	and return the trinity fq input string"""
	fq1,fq2 = [],[] #list of forward and reverse read files
	assert len(inDIRs) >= 1, "Empty input fasq.gz directories"
	print inDIRs
	for inDIR in inDIRs:
		if inDIR[-1] != "/": inDIR += "/"
		for i in os.listdir(inDIR):
			if inDIR == "./": inDIR = "" # So that the combined absolute dir look right
			if i.endswith(infile_end):
				if read1_identifier in i and read2_identifier not in i:
					fq1.append(inDIR+i)
					print "Adding forward read file",i
				elif read2_identifier in i and read1_identifier not in i:
					fq2.append(inDIR+i)
					print "Adding reverse read file",i
	assert len(fq1) > 0, "No forward read file found"
	assert len(fq2) > 0, "No reverse read file found"
	assert len(fq1) == len(fq2), "Unequal number of forward and reverse files"
	return " --left "+",".join(fq1)+" --right "+",".join(fq2)


def get_input_string_se(inDIRs,infile_end):
	"""Collect all the single end fastq files 
	and return the trinity fq input string"""
	fq = [] #list of single end read files
	assert len(inDIRs) >= 1, "Empty input fasq.gz directories"
	for inDIR in inDIRs: 
		if inDIR[-1] != "/": inDIR += "/"
		for i in os.listdir(inDIR):
			if i.endswith(infile_end):
				if inDIR == "./": inDIR = "" 
				fq.append(inDIR+i)
	for i in fq: print "Adding read file",i
	assert len(fq) > 0, "No fastq file found"
	return " --single "+",".join(fq)

	
def run_trinity(input_string,taxonID,num_cores,max_memory_GB,stranded=False,clip=True):
	"""Assemble using trinity v2"""
#	if clip:
#		assert os.path.exists(TruSeq_ADAPTER),"Cannot fine the adapter file "+TruSeq_ADAPTER
#		trim_setting = '"ILLUMINACLIP:'+TruSeq_ADAPTER+':2:30:10 SLIDINGWINDOW:4:5 LEADING:5 TRAILING:5 MINLEN:25"'
#	else: trim_setting = '"SLIDINGWINDOW:4:5 LEADING:5 TRAILING:5 MINLEN:25"'
	logfile=taxonID+".log"
	transcripts = taxonID+".Trinity.fasta"
	if stranded:
		strand = " --SS_lib_type RF"
	else: strand = ""
	if fasta_ok(transcripts) or os.path.exists(transcripts+".gz"):
		print "Skip trinity"
	else:
		cmd ="ulimit -s unlimited\n"
		cmd += TRINITY_CMD+" --seqType fq"
		#cmd += " --trimmomatic --quality_trimming_params "+trim_setting
		cmd += " --max_memory "+max_memory_GB+"G --CPU "+num_cores
		cmd += " --full_cleanup"
		cmd += " --no_normalize_reads"
		cmd += strand+" --output "+taxonID+".trinity"
		cmd += input_string
		os.system(TRINITY_CMD+" --version >>"+logfile) # log the version
		run(cmd,logfile)
		assert fasta_ok(taxonID+".trinity.Trinity.fasta"), \
			"Trinity did not finish correctly"
		os.rename(taxonID+".trinity.Trinity.fasta",transcripts) # shorten name


def shorten_fasta_names(inname,outname,taxonID):
	"""shorten transdecoder cds and pep file names"""
	infile = open(inname,"r")
	outfile = open(outname,"w")
	for line in infile:
		if line[0] == ">":
			newid = (line.split(" ")[0]).split(".")[-1]
			outfile.write(">"+taxonID+"@"+newid+"\n")
		else: outfile.write(line)
	infile.close()
	outfile.close()

				
def run_transdecoder_blastp(taxonID,num_cores,stranded=False,pfam=False):
	"""translate"""
	logfile=taxonID+".log"
	transcripts = taxonID+".Trinity.fasta"
	blastpout = taxonID+".blastp.outfmt6"
	if fasta_ok(taxonID+".pep.fa") and fasta_ok(taxonID+".cds.fa") and \
		os.path.exists(blastpout+".gz"):
		print "Skip transdecoder"
		return
	outpep,outcds = transcripts+".transdecoder.pep",transcripts+".transdecoder.cds"
	if fasta_ok(outpep) and fasta_ok(outcds): print "Skip transdecoder"
	else:		
		# Get all the candidate ORFs on the plus strand
		allpep = transcripts+".transdecoder_dir/longest_orfs.pep"
		if fasta_ok(allpep): print "Skip looking for long orfs"
		else:
			if stranded: strand = " -S"
			else: strand = ""
			cmd = TRANSDECODER_PFEFIX+"LongOrfs -t "+transcripts+strand
			run(cmd,logfile)
			assert fasta_ok(allpep), allpep+"small in size"
		
		# Compare to know peptides.
		# I use a blastp data set of Beta vulgaris and Arabidopsis thaliana
		if blastpout_ok(blastpout): print "Skip blastp"
		else:
			assert os.path.exists(BLASTP_DB_PATH),"cannot find "+BLASTP_DB_PATH
			cmd = "blastp -query "+allpep+" -db "+BLASTP_DB_PATH
			cmd += " -max_target_seqs 1 -outfmt 6 -evalue 10 -num_threads "+str(num_cores)
			cmd += " > "+blastpout
			run(cmd,logfile)
			assert blastpout_ok, "Few blastp hits"
			
		if pfam: # Compare to pfam A and B
			hmmout = taxonID+".pfam.domtblout"
			cmd = HMMSCAN_PATH+" --cpu "+str(num_cores)+" --domtblout "+hmmout+" "
			cmd += PFAM_PATH+" "+allpep
			run(cmd,logfile)
			with open(taxonID+".transdecoder-log","a") as outfile:
				outfile.write(cmd+"\n")
			
		# Get final peptides and CDS
		if fasta_ok(outpep) and fasta_ok(outcds): print "Skip finding final cds and pep"
		else:	
			cmd = TRANSDECODER_PFEFIX+"Predict -t "+transcripts
			if pfam: cmd += " --retain_pfam_hits "+hmmout
			cmd += " --retain_blastp_hits "+blastpout
			cmd += " --cpu "+str(num_cores)
			run(cmd,logfile)
		assert fasta_ok(outpep,min_count=500) and fasta_ok(outcds,min_count=500), \
			"transdecoder did not finish correctly"
	
	# shorten names
	shorten_fasta_names(outpep,taxonID+".pep.fa",taxonID)
	shorten_fasta_names(outcds,taxonID+".cds.fa",taxonID)
	# compress original output
	os.system("gzip "+outcds)
	os.system("gzip "+outpep)
	os.system("gzip "+blastpout)
	os.system("gzip "+transcripts)
	if pfam: os.system("gzip "+clusterID+"pfam.domtblout")
	
	try: # remove intermediate files
		os.remove(transcripts+".transdecoder.gff3")
		os.remove(transcripts+".transdecoder.mRNA")
		os.remove(transcripts+".transdecoder.bed")
		os.remove(transcripts+".transdecoder.dat")
		shutil.rmtree()(+transcripts+".transdecoder_dir")
	except: pass # ok if the intermediate files are removed already
	print "outfiles written to",taxonID+".pep.fa",taxonID+".cds.fa"


def check_pep_coverage_redundancy(taxonID,hitID="Beta",min_pident=60.0,log=True):
	"""blastp output tabular:
	0-qseqid	1-sseqid	2-pident	3-length	4-mismatch 	5-gapopen 
	6-qstart	7-qend		8-sstart 	9-send		10-evalue	11-bitscore
	Only summerize over hits with id starting with hitID, i.e. the closest reference
	"""
	blastpout = taxonID+".blastp.outfmt6.gz" # first column looks like TR7|c0_g1_i1|m.1
	cds = taxonID+".Trinity.fasta.transdecoder.pep.gz" # seq id looks like
	# >TR10002|c0_g1_i1|m.6695 TR10002|c0_g1_i1|g.6695  ORF TR10002|c0_g1_i1|g.6695 TR10002|c0_g1_i1|m.6695 type:complete len:104 (-) TR10002|c0_g1_i1:834-1145(-)
	logfile = taxonID+".log"
	cdsids = [] # all the ids in the cds file
	infile = gzip.open(cds,"rb")
	for line in infile:
		if line[0] == ">": cdsids.append((line.split(" ")[0])[1:])
	infile.close()
	assert len(cdsids) > 0, cds+"is not in correct fasta format"
	infile = gzip.open(blastpout,"rb")
	query_list = [] # list of all query ids
	hitcov_dict = {} # key is hit seqid, value is the max hit coverage
	for line in infile:
		spls = line.strip().split("\t")
		query,hit,piden = spls[0],spls[1],float(spls[2])
		if query not in cdsids:
			continue # only the ones get translated matters
			# Turns out that the cds file contains all seqs that has blastp hits
			# Filtering for ids in cds file doesn't really matter
		query_list.append(query)
		if piden < min_pident or not hit.startswith(hitID):
			continue # only look at highly similar hits to the closest proteome
		hitcov = int(spls[9])-int(spls[8])
		if hit not in hitcov_dict:
			hitcov_dict[hit] = hitcov
		else: hitcov_dict[hit] = max(hitcov_dict[hit],hitcov)
	infile.close()
	sum_hitcov_bp = 0
	for hitid in hitcov_dict:
		sum_hitcov_bp += hitcov_dict[hitid]
	sum_hitcov_numgen = len(hitcov_dict)
	out = "total reference coverage in bp: "+str(sum_hitcov_bp)+"\n"
	out += "total number of reference seqids covered: "+str(sum_hitcov_numgen)+"\n"
	out += "redundancy: "+str(len(set(query_list))/float(sum_hitcov_numgen))+"\n"
	if log:
		with open(logfile,"a") as outfile:
			outfile.write(out)
	print out
