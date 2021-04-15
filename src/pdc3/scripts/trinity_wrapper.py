"""
Runs Trinity
Triniy and Transrate must be in the path

"""

import os,sys
import shutil

def run(cmd,logfile):
	"""print, run and log calls"""
	print(cmd)
	os.system(cmd)
	with open(logfile,"a") as outfile: outfile.write(cmd+"\n")


def run_trinity_pe(pe_fq1,pe_fq2,taxonID,num_cores,max_memory_GB,strand,DIR):

	if DIR == ".": DIR = os.getcwd()	
	if os.path.isabs(DIR) == False: DIR = os.path.abspath(DIR)
	if DIR[-1] != "/": DIR += "/"

	logfile= taxonID+".log"
	transcripts = taxonID+".Trinity.fasta"
		
	if strand == "stranded":
		stranded = " --SS_lib_type RF"
	else: 
		stranded = ""

	if os.path.exists(DIR+transcripts):
		print("Transcript file found, skipping Trinity")
	else:
		cmd ="ulimit -s unlimited\n"
		cmd += "Trinity"+" --seqType fq"
		cmd += " --max_memory "+str(max_memory_GB)+"G --CPU "+str(num_cores)
		cmd += " --bflyCalculateCPU"
		cmd += " --full_cleanup"
		cmd += " --no_normalize_reads"
		cmd += stranded
		cmd += " --left "+pe_fq1
		cmd += " --right "+pe_fq2
		cmd += " --output "+DIR+taxonID+".trinity"
		os.system("Trinity"+" --version >>"+(DIR+logfile)) # log the version
		run(cmd,(DIR+logfile))
		assert (DIR+taxonID+".trinity.Trinity.fasta"), \
			"Trinity did not finish correctly"
		os.rename(DIR+taxonID+".trinity.Trinity.fasta",DIR+transcripts) # shorten name

def run_trinity_se(se_fq,taxonID,num_cores,max_memory_GB,strand,DIR):

	if DIR == ".": DIR = os.getcwd()	
	if os.path.isabs(DIR) == False: DIR = os.path.abspath(DIR)
	if DIR[-1] != "/": DIR += "/"

	logfile= taxonID+".log"
	transcripts = taxonID+".Trinity.fasta"
		
	if strand == "stranded":
		stranded = " --SS_lib_type R"
	else: 
		stranded = ""

	if os.path.exists(DIR+transcripts):
		print("Transcript file found, skipping Trinity")
	else:
		cmd ="ulimit -s unlimited\n"
		cmd += "Trinity"+" --seqType fq"
		cmd += " --max_memory "+str(max_memory_GB)+"G --CPU "+str(num_cores)
		cmd += " --bflyCalculateCPU"
		cmd += " --full_cleanup"
		cmd += " --no_normalize_reads"
		cmd += stranded
		cmd += " --single "+se_fq
		cmd += " --output "+DIR+taxonID+".trinity"
		os.system("Trinity"+" --version >>"+(DIR+logfile)) # log the version
		run(cmd,(DIR+logfile))
		assert (DIR+taxonID+".trinity.Trinity.fasta"), \
			"Trinity did not finish correctly"
		os.rename(DIR+taxonID+".trinity.Trinity.fasta",DIR+transcripts) # shorten name
		
if __name__ == "__main__":

	if len(sys.argv) == 8:
		pe_fq1,pe_fq2,taxonID,num_cores,max_memory_GB,strand,DIR = sys.argv[1:]
		assert strand == "stranded" or strand == "non-stranded", \
		"strand has to be either stranded or non-stranded"
		run_trinity_pe(pe_fq1,pe_fq2,taxonID,num_cores,max_memory_GB,strand,DIR)
		
	elif len(sys.argv) == 7:
		se_fq,taxonID,num_cores,max_memory_GB,strand,DIR = sys.argv[1:]
		assert strand == "stranded" or strand == "non-stranded", \
		"strand has to be either stranded or non-stranded"
		run_trinity_se(se_fq,taxonID,num_cores,max_memory_GB,strand,DIR)	
	else:
		print ("Usage:")
		print ("For single end reads: python trinity_wrapper.py fastq_se_reads taxonID num_cores max_memory_GB strand output_dir")
		print ("For paired end reads: python trinity_wrapper.py fastq_pe_reads1 fastq_pe_reads2 taxonID num_cores max_memory_GB strand(stranded or non-stranded) output_dir")
		sys.exit(0)


