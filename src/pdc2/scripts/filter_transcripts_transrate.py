"""
Filter low quality transcripts indenfided by Transrate.
It uses individual cutoffs for sCord, sCcov and
sCnuc rather that transrate contig scores.

"""

import sys, os, re
import pandas as pd
from Bio import SeqIO


def filter_bad_transcripts(transcripts,transrate_contings_csv,DIR):

	if DIR == ".": DIR = os.getcwd()
	if os.path.isabs(DIR) == False: DIR = os.path.abspath(DIR)
	if DIR[-1] != "/": DIR += "/"

	path_transcript, files_transcript = os.path.split(transcripts)
	transcripts_name = str(files_transcript)
	base_name_transcripts = transcripts_name.split( "." )
	
	good_transcripts_file = base_name_transcripts[0]+".good_transcripts.fa"
	bad_transcripts_file = base_name_transcripts[0]+".bad_transcripts.fa"
	good_transcripts_short_name_file = base_name_transcripts[0]+".good_transcripts.short_name.fa"
    
	if os.path.exists(DIR+good_transcripts_file) and os.path.exists(DIR+bad_transcripts_file) \
		and os.path.exists(DIR+good_transcripts_short_name_file):
		print "Filter transcript files found for", base_name_transcripts[0]
	
	else:
    
    	#open transrate_contig_csv and transcript fasta files
		contigs = pd.read_csv(transrate_contings_csv)
		transcripts_original = SeqIO.index(transcripts, "fasta")
	
	
    	# define cutoff values, change if needed
		Cord = contigs["sCord"]
		Cord_cutoff = Cord<=0.50
	
		Ccov = contigs["sCcov"]
		Ccov_cutoff = Ccov<=0.25
	
		Cnuc = contigs["sCnuc"]
		Cnuc_cutoff = Cnuc<=0.25
	
    	#filer transrate_contig_csv by component
		misassembly = contigs[Cord_cutoff]
		uncovered = contigs[Ccov_cutoff]
		nonagreement = contigs[Cnuc_cutoff]
	
    	#merge all bad transcripts, remove duplicated row and get sequences names
		components = [misassembly,uncovered,nonagreement]
		reads_to_filter = pd.concat(components)
		undup_reads_to_filter = pd.DataFrame.drop_duplicates(reads_to_filter, keep='first')
		bad_transcripts_names = undup_reads_to_filter["contig_name"]
	
    	#write fasta and csv files of bad transcripts
		bad_transcripts= []
		if len(bad_transcripts_names) == 0:
			print "Bad transcripts not found"
		for i in bad_transcripts_names: bad_transcripts.append(transcripts_original[i])
		count = SeqIO.write(bad_transcripts,(DIR+bad_transcripts_file), "fasta")
		print("Removed %i bad transcripts" % count)
		undup_reads_to_filter.to_csv(DIR+str(base_name_transcripts[0])+".bad_transcripts.csv",index=False)
	
    	#filter only good transcripts and sequences names
	
		good_transcripts_to_keep = pd.concat([contigs,undup_reads_to_filter]).drop_duplicates(keep=False)
		good_transcripts_names = good_transcripts_to_keep["contig_name"]
		good_transcripts = []
	
    	#write fasta and csv files of good transcripts
		if len(good_transcripts_names) == 0:
			print "Good transcripts not found"
		for i in good_transcripts_names: good_transcripts.append(transcripts_original[i])
		count = SeqIO.write(good_transcripts,(DIR+good_transcripts_file), "fasta")
		print("Kept %i good transcripts" % count)
		good_transcripts_to_keep.to_csv(DIR+str(base_name_transcripts[0])+".good_transcripts.csv",index=False)
	
	
		#trim trinity length and path information from transcript fasta names of good contigs
		searchstr = r'(>\w+)(\slen.*)'
		replacestr = r'\1'
		infile = open((DIR+good_transcripts_short_name_file), 'w')
	
		with open((DIR+good_transcripts_file), 'rU') as fasta_file:
			reg = re.compile(searchstr)
			for line in fasta_file:
				line=line.strip('\n')
				if line.startswith('>'):
					fixline = reg.sub(replacestr, line)
					infile.write(fixline + '\n')
				else:
					infile.write(line + '\n')
	
		infile.close()
		
	assert os.path.exists(DIR+good_transcripts_file) and os.path.exists(DIR+bad_transcripts_file) \
		and os.path.exists(DIR+good_transcripts_short_name_file), "filter_transcript_transrate did not finish"


if __name__ == "__main__":
	if len(sys.argv) == 4:
		filter_bad_transcripts(transcripts=sys.argv[1],transrate_contings_csv=sys.argv[2],DIR=sys.argv[3])
	else:
		print ("Usage:")
		print ("python filter_transcripts_transrate.py transcripts_fasta transrate_csv output_dir")
		sys.exit(0)
