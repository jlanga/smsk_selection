"""
Wrapper for Rcorrector. 
It takes either single or paired-end fastq files
Rcorrector needs to be installed and in the path
The output are fastq files with fixed reads flagged as 'cor'
and unfixable reads flagged as 'unfixable_error'
It can handle .gz files
"""

import sys,os

APPS_HOME = "/usr/local/bin/" # where rcorrector is located
RCORRECTOR_CMD = 'perl '+APPS_HOME+"/Rcorrector/run_rcorrector.pl" # Basic command to call rcorrector


def rcorrector_se(se_fq,num_cores,DIR):

	if DIR == ".": DIR = os.getcwd()	
	if os.path.isabs(DIR) == False: DIR = os.path.abspath(DIR)
	if DIR[-1] != "/": DIR += "/"

	path_se, file_se = (os.path.split(se_fq)) #splits the path from the file name
	se_name = str(file_se)
	base_name_se = se_name.split( "." )
	
	if (os.path.splitext(file_se)[-1]) == ".gz" :	
		corrected = (base_name_se[0])+".cor.fq.gz"
	else:
		corrected = (base_name_se[0])+".cor.fq"
	
	if os.path.exists(DIR+corrected): 
		print ("Found", corrected)
	else:
		cmd = [RCORRECTOR_CMD,"-s",se_fq,"-t",str(num_cores),"-od",DIR]
		print (" ".join(cmd))
		os.system(" ".join(cmd))
	assert os.path.exists(DIR+corrected),"Rcorrector did not finished"

def rcorrector_pe(pe_fq1,pe_fq2,num_cores,DIR):

	if DIR == ".": DIR = os.getcwd()	
	if os.path.isabs(DIR) == False: DIR = os.path.abspath(DIR)
	if DIR[-1] != "/": DIR += "/"

	path_pe_1, file_pe_1 = (os.path.split(pe_fq1)) #splits the path from the file name
	path_pe_2, file_pe_2 = (os.path.split(pe_fq2))

	path_pe_1, file_pe_1 = os.path.split(pe_fq1)
	pe_fq1_name = str(file_pe_1)
	base_name_1 = pe_fq1_name.split( "." )
	
	path_pe_2, file_pe_2 = os.path.split(pe_fq2)
	pe_fq2_name = str(file_pe_2)
	base_name_2 = pe_fq2_name.split( "." )

			
	if (os.path.splitext(file_pe_1)[-1]) and (os.path.splitext(file_pe_2)[-1]) == ".gz" :	
		corrected_1 = (base_name_1[0])+".cor.fq.gz" 
		corrected_2 = (base_name_2[0])+".cor.fq.gz" 
	else:
		corrected_1 = (base_name_1[0])+".cor.fq" #input paired fq1 with extension ".cor.fq"
		corrected_2 = (base_name_2[0])+".cor.fq" #input paired fq2 with extension ".cor.fq"
	
	
	if os.path.exists(DIR+corrected_1) and os.path.exists(DIR+corrected_2): 
		print ("Found", corrected_1, corrected_2)
	else:
		cmd = [RCORRECTOR_CMD,"-1",pe_fq1,"-2",pe_fq2,"-t",str(num_cores),"-od",DIR]
		print (" ".join(cmd))
		os.system(" ".join(cmd))
		
	assert os.path.exists(DIR+corrected_1) \
			and os.path.exists(DIR+corrected_2),"Rcorrector did not finished"

if __name__ == "__main__":
	if len(sys.argv) == 4:
		rcorrector_se(se_fq=sys.argv[1],num_cores=int(sys.argv[2]),DIR=sys.argv[3])
	elif len(sys.argv) == 5:
		rcorrector_pe(pe_fq1=sys.argv[1],pe_fq2=sys.argv[2],num_cores=int(sys.argv[3]),DIR=sys.argv[4])
	else:
		print ("Usage:")
		print ("For single end reads: python rcorrector.py fastq_se_reads num_cores output_dir")
		print ("For paired end reads: python rcorrector.py fastq_pe_reads1 num_cores fastq_pe_reads2 output_dir")
		sys.exit(0)
		
		


