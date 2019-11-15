"""
All reads generated in Smith Lab are stranded and paired-end 101 bp read from HiSeq.
Read trimming is carried out using trimmomatic with custome adaptor files. 
"""

import os,sys
import process_seq

if __name__ =="__main__":
	if len(sys.argv) != 5:
		print("usage: python run_all_illumina_pe_stranded.py fastq.gzDIR taxonID num_cores max_memory_GB")
		print("cd into the output directory.")
		print("separate by ',' where there are multiple fastq.gz directories")
		sys.exit()
	
	DIR,taxonID,num_cores,max_memory_GB = sys.argv[1:5]
	
	if not os.path.exists(taxonID+".Trinity.fasta"):
		# The original fastq reads downloaded form the sequqncing core looks like
		# 41620_TGACCA_L005_R1_001.fastq.gz
		# 41620_TGACCA_L005_R2_001.fastq.gz
		# With one directory per sample
		# Collect all the fastq.gz files and return the trinity fq input string
		fqstring = process_seq.get_input_string_pe(DIR.split(","),\
			read1_identifier="_R1_",\
			read2_identifier="_R2_",\
			infile_end=".fastq.gz")
	
		# assemble
		process_seq.run_trinity(fqstring,taxonID,num_cores,max_memory_GB,\
			stranded=True,clip=True)
	
	# translate
	process_seq.run_transdecoder_blastp(taxonID,num_cores,stranded=True)
	
	# summarize coverage and redundancy
	process_seq.check_pep_coverage_redundancy(taxonID)
