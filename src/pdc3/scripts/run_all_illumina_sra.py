"""
fastq-dump all .sra files, combine and assemble as one data set

Given the blackbox nature and change of phred scores
Use the olc filtering scirpts

Change FASTQ_DUMP_CMD and GENERAL_ADAPTER as needed
"""

import os,sys
import process_seq
import filter_fastq

FASTQ_DUMP_CMD = "fastq-dump.2.8.2" # part of the sra toolkit

# GENERAL_ADAPTER contains both TruSeq plus various other adapters
# Use this instead of TruSeq_ADAPTER if using SRA or reads of uncertain sources
GENERAL_ADAPTER = "/home/yayang/phylogenomic_dataset_construction/data/UniVec-TruSeq_adapters"


if __name__ =="__main__":
	if len(sys.argv) != 7:
		print("usage: python run_all_illumina_sra.py .sraDIR taxonID single/pair stranded/non-stranded num_cores max_memory_GB")
		print("each .sra DIR will be assembled as one data set")
		sys.exit()
	
	DIR, taxonID, seqtype, strand, num_cores, max_memory_GB = sys.argv[1:]
	assert seqtype == "single" or seqtype == "pair", \
		"sequence type has to be either single or pair"
	assert strand == "stranded" or strand == "non-stranded", \
		"strand has to be either stranded ornon-stranded"
	logfile = taxonID+".log"
	if DIR[-1] != "/": DIR += "/"
	
	# extract sra archives to fastq
	sra_names = [] # list of sra file names minus .sra
	for i in os.listdir(DIR):
		if i.endswith(".sra"): # file looks like SRR1559276.sra
			sra_names.append(DIR+i.replace(".sra",""))
	assert len(sra_names) > 0, "No file end with .sra found in "+DIR
	for name in sra_names:
		cmd = "fastq-dump.2.8.2 --defline-seq '@$sn[_$rn]/$ri' "+name+".sra --outdir "+DIR[:-1]
		if seqtype == "pair":
			rawfq1,rawfq2 = name+"_1.fastq",name+"_2.fastq"
			if filter_fastq.fastq_ok(rawfq1) and filter_fastq.fastq_ok(rawfq2):
				print("Found", rawfq1, rawfq2)
			else:
				cmd += " --split-files"
				process_seq.run(cmd,logfile)
				assert filter_fastq.fastq_ok(rawfq1) and filter_fastq.fastq_ok(rawfq2),\
					"fastq-dump did not finish correctly"
			filter_fastq.filter_fastq_pe(rawfq1,rawfq2,GENERAL_ADAPTER,num_cores,logfile)
		else:
			rawfq = name+".fastq"
			if filter_fastq.fastq_ok(rawfq):
				print("Found", rawfq)
			else:
				process_seq.run(cmd,logfile)
				assert filter_fastq.fastq_ok(rawfq),"fastq-dump did not finish correctly"
			filter_fastq.filter_fastq_se(rawfq,GENERAL_ADAPTER,num_cores,logfile)
		
	#assemble
	stranded = True if strand == "stranded" else False
	if seqtype == "pair":
		fqstring = process_seq.get_input_string_pe(inDIRs=[DIR],\
				   read1_identifier="_1.",\
				   read2_identifier="_2.",\
				   infile_end=".fastq.filtered")
	else: fqstring = process_seq.get_input_string_se(inDIRs=[DIR],\
				   infile_end=".fastq.filtered")
	process_seq.run_trinity(fqstring,taxonID,num_cores,max_memory_GB,stranded=stranded,clip=False)
	
	# translate
	process_seq.run_transdecoder_blastp(taxonID,num_cores,stranded=stranded)
	
	# summarize coverage and redundancy
	process_seq.check_pep_coverage_redundancy(taxonID)
	
	os.system("gzip "+DIR+"*remove")

