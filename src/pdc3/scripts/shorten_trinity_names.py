"""
trim trinity length and path information from transcript fasta names of good contigs
"""

import sys, os, re

def short_trinity_names(transcripts,DIR):

	if DIR == ".": DIR = os.getcwd()
	if os.path.isabs(DIR) == False: DIR = os.path.abspath(DIR)
	if DIR[-1] != "/": DIR += "/"

	path_transcript, files_transcript = os.path.split(transcripts)
	transcripts_name = str(files_transcript)
	base_name_transcripts = transcripts_name.split( "." )

	transcripts_short_name_file = base_name_transcripts[0]+".Trinity.short_name.fa"

	if DIR == ".": DIR = os.getcwd()
	if os.path.isabs(DIR) == False: DIR = os.path.abspath(DIR)
	if DIR[-1] != "/": DIR += "/"

	#trim trinity length and path information from transcript fasta names of good contigs
	searchstr = r'(>\w+)(\slen.*)'
	replacestr = r'\1'
	outfile = open((DIR+transcripts_short_name_file), 'w')
	
	with open((DIR+transcripts), 'rU') as fasta_file:
		reg = re.compile(searchstr)
		for line in fasta_file:
			line=line.strip('\n')
			if line.startswith('>'):
				fixline = reg.sub(replacestr, line)
				outfile.write(fixline + '\n')
			else:
				outfile.write(line + '\n')
	
	outfile.close()
	
	
if __name__ == "__main__":
	if len(sys.argv) == 3:
		short_trinity_names(transcripts=sys.argv[1],DIR=sys.argv[2])
	else:
		print ("Usage:")
		print ("python shorten_trinity_names.py transcripts_fasta output_dir")
		sys.exit(0)

