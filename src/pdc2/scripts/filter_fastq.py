"""
Compare paired-end illumina reads with user-supplied adapter sequence
If adapters are detected in either reads, remove the entire read pair
Filter reads by quality scores before output .fq.filtered files

Dependancies: makeblastdb, blastn
"""

import sys,os
import seq

SAMPLE_FREQ = 1000000

def run_cmd(cmd,logfile):
	"""print, log and run"""
	print cmd
	with open(logfile,"a") as infile:
		infile.write(cmd+"\n")
	os.system(cmd)

def fastq_ok(fqfile):
	"""Check a fastq file is in the right format"""
	if not os.path.exists(fqfile): return False
	count = 0
	infile = open(fqfile,"r")
	for read_obj in seq.fastq_generator(infile):
		count += 1
		if count == 100: break
	infile.close()
	if count == 100: return True
	else: return False


def blast_adapters(fq,adapter_file,num_cores,logfile):
	"""identify reads with adapter contamination"""
	fa = fq+".fasta"
	if os.path.exists(fq+".to_remove") and os.stat(fq+".to_remove").st_size > 0:
		return fq+".to_remove"
	#convert fastq to fasta
	if not os.path.exists(fa):
		cmd = "sed -n '1~4s/^@/>/p;2~4p' "+fq+" > "+fa
		run_cmd(cmd,logfile)
	#blast and pipe the output hits
	os.system("makeblastdb -in "+fa+" -parse_seqids -dbtype nucl -out "+fa+".db")
	cmd = "blastn -query "+adapter_file
	cmd += " -db "+fa+".db -outfmt 6 -num_threads "+str(num_cores)+" -max_target_seqs 100000000 "
	cmd += "| awk '{print $2}' | sort | uniq >"+fq+".to_remove"
	run_cmd(cmd,logfile)
	if os.path.exists(fq+".to_remove") and os.stat(fq+".to_remove").st_size > 0:
		os.system("rm "+fa+".db* "+fa)
	else: sys.exit("Error in blastn")
	return fq+".to_remove"


def filter_fastq_se(fq,adapter_file,num_cores,logfile="log"):
	"""Filter single end reads"""
	out = fq+".filtered"
	if fastq_ok(out):
		print out,"exists"
		return out
	#detect adapters
	to_remove = blast_adapters(fq,adapter_file,num_cores,logfile)
	with open(to_remove) as infile:
		adapterid_set = set(infile.read().splitlines())
	print len(adapterid_set),"reads will be removed"
	before, after = 0, 0
	infile = open(fq,"r")
	outfile = open(out,"w")
	outfile_qual = open(fq+".qual","w") # quality scores for plotting in R
	for read_obj in seq.fastq_generator(infile):
		before += 1 #keep track of number of input read pairs
		#shorten names here to match blast hits processed by awk
		name = (read_obj.name).split(" ")[0]
		if name in adapterid_set:
			#remove seqid to speed up remaining searches
			adapterid_set.remove(name)
		else: #does not contain adapter
			after += 1
			outfile.write(read_obj.get_fastq())		
		if before % SAMPLE_FREQ == 0:
			outfile_qual.write(" ".join([str(score) for score in read_obj.qualarr])+"\n")
			print "Read",before,"reads,",after,"written to outfiles"
	outfile.close()
	outfile_qual.close()
	infile.close()
	with open(logfile,"a") as outfile:
		outfile.write("input file: "+fq+"\n")
		outfile.write("adapter file: "+adapter_file+"\n")
		outfile.write("Original reads: "+str(before)+"\n")
		outfile.write("Reads after removing adapters: "+str(after)+"\n")
	print "Output written to .filtered files"
	print "Summary stats written to",logfile
	assert fastq_ok(out), "Error in filtering fastq files"
	return out


def filter_fastq_pe(fq1,fq2,adapter_file,num_cores,logfile="log"):
	"""Filter paired end reads"""
	out1,out2 = fq1+".filtered",fq2+".filtered"
	if fastq_ok(out1) and fastq_ok(out2):
		print out1,out2,"exists"
		return out1,out2
	#detect adapters
	to_remove1 = blast_adapters(fq1,adapter_file,num_cores,logfile)
	to_remove2 = blast_adapters(fq2,adapter_file,num_cores,logfile)
	with open(to_remove1) as infile:
		ids = set(infile.read().splitlines())
	with open(to_remove2) as infile:
		ids_to_remove = ids.union(set(infile.read().splitlines()))
	print len(ids_to_remove),"read pairs will be removed"
	print "Writing the processed fastq files"
	
	before, after = 0, 0
	infile = open(fq1,"r")
	outfile = open(out1,"w")
	outfile_qual = open(fq1+".qual","w") # quality scores for plotting in R
	for read_obj in seq.fastq_generator(infile):
		before += 1
		# cut names here to match blast hits processed by awk
		name = (read_obj.name).split(" ")[0]
		if name not in ids_to_remove:
			# add /1 and /2 to be compatible to Trinity
			read_obj.name = name+"/1"
			after += 1
			outfile.write(read_obj.get_fastq())
		if before % SAMPLE_FREQ == 0:
			outfile_qual.write(" ".join([str(score) for score in read_obj.qualarr])+"\n")
			print "Read",before,"pairs,",after,"written to outfiles"
	outfile.close()
	outfile_qual.close()
	infile.close()
	with open(logfile,"a") as outfile:
		outfile.write("adapter file: "+adapter_file+"\n")
		outfile.write("Original reads: "+str(before)+"\n")
		outfile.write("Reads after removing adapters: "+str(after)+"\n")
	
	#reset for the reverse reads
	before, after = 0, 0
	infile = open(fq2,"r")
	outfile = open(out2,"w")
	outfile_qual = open(fq2+".qual","w") # quality scores for plotting in R
	for read_obj in seq.fastq_generator(infile):
		before += 1
		# cut names here to match blast hits processed by awk
		name = (read_obj.name).split(" ")[0]
		if name not in ids_to_remove:
			# add /1 and /2 to be compatible to Trinity
			read_obj.name = name+"/2"
			after += 1
			outfile.write(read_obj.get_fastq())
		if before % SAMPLE_FREQ == 0:
			outfile_qual.write(" ".join([str(score) for score in read_obj.qualarr])+"\n")
			print "Read",before,"pairs,",after,"written to outfiles"
	outfile.close()
	outfile_qual.close()
	infile.close()
	with open(logfile,"a") as outfile:
		outfile.write("adapter file: "+adapter_file+"\n")
		outfile.write("Original reads: "+str(before)+"\n")
		outfile.write("Reads after removing adapters: "+str(after)+"\n")		
		
	print "Output written to .filtered files"
	print "Summary stats written to",logfile
	assert fastq_ok(out1) and fastq_ok(out2), "Error in filtering fastq files"
	return out1,out2

if __name__ == "__main__":
	if len(sys.argv) == 4:
		filter_fastq_se(fq=sys.argv[1],adapter_file=sys.argv[2],num_cores=int(sys.argv[3]))
	elif len(sys.argv) == 5:
		filter_fastq_pe(fq1=sys.argv[1],fq2=sys.argv[2],adapter_file=sys.argv[3],num_cores=int(sys.argv[4]))
	else:
		print "Usage:"
		print "For single end reads: python filter_fastq.py fq_read adapter_file num_cores"
		print "For paired end reads: python filter_fastq.py fq_read1 fq_read2 adapter_file num_cores"
		sys.exit(0)
	


