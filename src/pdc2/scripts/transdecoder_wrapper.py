"""
Transdecoder wrapper
transdecoder and blastp must be in the path
"""

import os

MIN_FASTA = 1000 # miminal number of seqs expected in assembly and translation files
BLASTP_DB_PATH = os.path.expanduser("~/Desktop/botany_2018/databases/db") # the custom blast database
CDS_REFERENCE = os.path.expanduser("~/Desktop/botany_2018/databases/Beta.fa")


import os,sys,shutil
from shutil import copyfile
import seq
import gzip

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


def shorten_fasta_names(inname,outname,taxonID):
	"""shorten transdecoder cds and pep file names"""
	infile = open(inname,"r")
	outfile = open(outname,"w")
	for line in infile:
		if line[0] == ">":
			newid = (line.split(" ")[0]).split(".")[0].replace(">TRINITY_","")
			outfile.write(">"+taxonID+"@"+newid+"\n")
		else: outfile.write(line)
	infile.close()
	outfile.close()

				
def run_transdecoder_blastp(transcripts,num_cores,strand,DIR):
	"""translate"""

	if DIR == ".": DIR = os.getcwd()	
	if os.path.isabs(DIR) == False: DIR = os.path.abspath(DIR)
	if DIR[-1] != "/": DIR += "/"

	path_transcript, files_transcript = os.path.split(transcripts)
	transcripts_name = str(files_transcript)
	base_name_transcripts = transcripts_name.split( "." )
	taxonID = base_name_transcripts[0]

	logfile = taxonID+".log"
	blastpout = taxonID+".blastp.outfmt6"

	if os.path.exists(DIR+files_transcript) == False:
		copyfile((os.path.abspath(transcripts)), (DIR+files_transcript))
	if DIR != os.getcwd():
		os.chdir(DIR)


	if (DIR+taxonID+".pep.fa") and (DIR+taxonID+".cds.fa") and \
		os.path.exists(DIR+blastpout+".gz"):
		print "Skip transdecoder"
		return

	outpep,outcds = files_transcript+".transdecoder.pep",files_transcript+".transdecoder.cds"
	if fasta_ok(DIR+outpep) and fasta_ok(DIR+outcds): print "Skip transdecoder"
	
	
	else:		
		# Get all the candidate ORFs on the plus strand
		allpep = files_transcript+".transdecoder_dir/longest_orfs.pep"

		if os.path.exists(DIR+allpep): 
			print "Skip looking for long orfs"

		else:
			if strand == "stranded":
				stranded = " -S"
			else: 
				stranded = ""
			
			cmd = "TransDecoder.LongOrfs -t "+(DIR+files_transcript)+stranded
			run(cmd,(DIR+logfile))
			assert fasta_ok(DIR+allpep), DIR+allpep+" small in size"	
	
		
		# Compare to know peptides.
		# I use a blastp data set of Beta vulgaris and Arabidopsis thaliana
		if blastpout_ok(DIR+blastpout): print "Skip blastp"
		else:
			assert os.path.exists(BLASTP_DB_PATH),"cannot find "+BLASTP_DB_PATH+ "pep database"
			
			cmd = "blastp -query "+(DIR+allpep)+" -db "+BLASTP_DB_PATH
			cmd += " -max_target_seqs 1 -outfmt 6 -evalue 10 -num_threads "+str(num_cores)
			cmd += " > "+(DIR+blastpout)
			run(cmd,(DIR+logfile))
			assert blastpout_ok(DIR+blastpout), DIR+blastpout+" with few hits"
			
			
		# Get final peptides and CDS
		if fasta_ok(DIR+outpep) and fasta_ok(DIR+outcds): print "Skip finding final cds and pep"
		else:	
			cmd = "TransDecoder.Predict -t "+(DIR+files_transcript)
			cmd += " --retain_blastp_hits "+(DIR+blastpout)
			cmd += " --cpu "+str(num_cores)
			run(cmd,(DIR+logfile))
		assert fasta_ok((DIR+outpep),min_count=500) and fasta_ok((DIR+outcds),min_count=500), \
			"transdecoder did not finish correctly"
	
	# shorten names
	shorten_fasta_names((DIR+outpep),(DIR+taxonID+".pep.fa"),taxonID)
	shorten_fasta_names((DIR+outcds),(DIR+taxonID+".cds.fa"),taxonID)
	# compress original output
	os.system("gzip "+(DIR+outpep))
	os.system("gzip "+(DIR+outcds))
	os.system("gzip "+(DIR+blastpout))
	#os.system("gzip "+(DIR+files_transcript))


	try: # remove intermediate files
		os.remove(DIR+files_transcript+".transdecoder.gff3")
		os.remove(DIR+files_transcript+".transdecoder.mRNA")
		os.remove(DIR+files_transcript+".transdecoder.bed")
		os.remove(DIR+files_transcript+".transdecoder.dat")
		shutil.rmtree(DIR+files_transcript+".transdecoder_dir")
		shutil.rmtree(DIR+files_transcript+".transdecoder_dir.__checkpoints")
	except: pass # ok if the intermediate files are removed already
	print "outfiles written to",(DIR+taxonID+".pep.fa"),(DIR+taxonID+".cds.fa")


def check_pep_coverage_redundancy(blastpout,hitID="Beta",min_pident=60.0,log=True):
	"""blastp output tabular:
	0-qseqid	1-sseqid	2-pident	3-length	4-mismatch 	5-gapopen 
	6-qstart	7-qend		8-sstart 	9-send		10-evalue	11-bitscore
	Only summerize over hits with id starting with hitID, i.e. the closest reference
	"""

	path_transcript, files_transcript = os.path.split(transcripts)
	transcripts_name = str(files_transcript)
	base_name_transcripts = transcripts_name.split( "." )
	taxonID = base_name_transcripts[0]
	
	cds = files_transcript+".transdecoder.pep.gz" # seq id looks like
	
	# >TRINITY_DN999_c0_g1_i1.p1 TRINITY_DN999_c0_g1~~TRINITY_DN999_c0_g1_i1.p1  ORF type:internal len:484 (+),score=73.41,Beta@Bv6_131850_txhn_t1|79.45|0.0 TRINITY_DN999_c0_g1_i1:2-1450(+)
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

	
def transrate_cds_ref(outcds,num_cores,reference=CDS_REFERENCE):
	"""Using transrate to get basic stats and reference coverage from final cds"""
	
	logfile = taxonID+".log"
	os.system("gunzip "+outcds)
	path_cds, file_cds = (os.path.split(outcds)) #splits the path from the file name
	cds_base_name = (os.path.splitext(file_cds)[0])
	results_name = cds_base_name+"translated_Transrate_results"
	
	if os.path.exists(results_name):
	
		print "Found cds transrate results folder"
	
	else:
	
		cmd = "transrate "+"--assembly "+cds_base_name
		cmd += " --threads "+str(num_cores)
		cmd += " --reference "+reference
		cmd += " --output "+results_name
		cmd += " > "+cds_base_name+"_translated_Transrate.log" #separate log for transrate for easy access
		run(cmd,logfile)

	assert os.path.exists(results_name), "Transrate not completed"
	
	os.system("gzip "+cds_base_name)

	
if __name__ == "__main__":

	if len(sys.argv) == 5:
	
		transcripts,num_cores,strand,DIR = sys.argv[1:]
		
		assert strand == "stranded" or strand == "non-stranded", \
			"strand has to be either stranded or non-stranded"
	
		path_transcript, files_transcript = os.path.split(transcripts)
		transcripts_name = str(files_transcript)
		base_name_transcripts = transcripts_name.split( "." )
		taxonID = base_name_transcripts[0]
	
		blastpout = taxonID+".blastp.outfmt6.gz" # first column looks like TR7|c0_g1_i1|m.1
		outcds = files_transcript+".transdecoder.cds.gz"
		
		run_transdecoder_blastp(transcripts,num_cores,strand,DIR)
		check_pep_coverage_redundancy(blastpout)
		#transrate_cds_ref(outcds,num_cores)

	else:
		print ("Usage:")
		print ("python transdecoder_wrapper.py transcripts num_cores strand(stranded or non-stranded) output_dir")
		sys.exit(0)
	
	
	

