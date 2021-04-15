"""
Wrapper to run blastx and chimera detection using method from 
Yang, Y. and S.A. Smith Optimizing de novo assembly of 
short-read RNA-seq data for phylogenomics. BMC Genomics 2013

blastx should be in the path

The output is a fasta file with the chimeras complete removed and a files with only chimeras 
(If you want to cut and keep the chimera pieces you need to run the cut_chimera_from_blastx.py
from the paper above)

"""

import sys,os
import pandas as pd
from Bio import SeqIO

SCRIPTS_HOME = os.path.expanduser("~/Desktop/botany_2018/scripts/") # where all scripts are located

#This make a database for blast given a reference proteome
def make_blast_db(proteome_ref,DIR):


	if DIR == ".": DIR = os.getcwd()	
	if os.path.isabs(DIR) == False: DIR = os.path.abspath(DIR)
	if DIR[-1] != "/": DIR += "/"


	path_ref, file_ref = os.path.split(proteome_ref) #splits the path from the file name
	ref_name = str(file_ref)
	ref_base_name = ref_name.split( "." )

	db_file_1 = ref_base_name[0]+".phr"
	db_file_2 = ref_base_name[0]+".pin"
	db_file_3 = ref_base_name[0]+".psq"


	if os.path.exists(DIR+db_file_1) and os.path.exists(DIR+db_file_2) \
		and os.path.exists(DIR+db_file_3): \
		print("BLAST database files found for", ref_base_name[0])
		
	else:	
		cmd = ["makeblastdb","-in",proteome_ref,"-dbtype prot","-out",(DIR+ref_base_name[0])]
		print((" ".join(cmd)))
		os.system(" ".join(cmd))	
		
	assert os.path.exists(DIR+db_file_1) and os.path.exists(DIR+db_file_2) \
		and os.path.exists(DIR+db_file_3), "makeblastdb did not finish"
	

#This runs a custom output of blastx
def run_blastx(transcripts,blast_db,num_cores,DIR):

	if DIR == ".": DIR = os.getcwd()	
	if os.path.isabs(DIR) == False: DIR = os.path.abspath(DIR)
	if DIR[-1] != "/": DIR += "/"

	path_transcript, file_transcript = os.path.split(transcripts) #splits the path from the file name
	transcript_name = str(file_transcript)
	transcript_base_name = transcript_name.split( "." )


	blastx_output_file = transcript_base_name[0]+".blastx"
	
	if os.path.exists(DIR+blastx_output_file):
		print("blastx output file found for",transcript_base_name[0])
	else:
		cmd = ["blastx","-db",blast_db,"-query",transcripts, \
		"-evalue 0.01","-outfmt",'"6 qseqid qlen sseqid slen frames pident nident length mismatch gapopen qstart qend sstart send evalue bitscore"', \
		"-out",(DIR+blastx_output_file),"-num_threads",str(num_cores),"-max_target_seqs 100"]
		print((" ".join(cmd)))
		os.system(" ".join(cmd))
			
	assert os.path.exists(DIR+blastx_output_file), "blastx did not finish"
	
	
#Call for detect_chimera_from_blastx_modifed.py
def run_chimera_detection(blastx_output_file,DIR):

	if DIR == ".": DIR = os.getcwd()	
	if os.path.isabs(DIR) == False: DIR = os.path.abspath(DIR)
	if DIR[-1] != "/": DIR += "/"

	path_blastx, file_blastx = os.path.split(blastx_output_file) #splits the path from the file name
	blastx_name = str(file_blastx)
	blastx_base_name = blastx_name.split( "." )

	cut_file = (blastx_base_name[0]+".cut")
	info_file = (blastx_base_name[0]+".info")

	if os.path.exists(DIR+cut_file) and os.path.exists(DIR+info_file):
		print("Chimera detection files found for", blastx_base_name[0])
		
	else:
		#detect_chimera_from_blastx_modifed(blastx_output_file,DIR)
		cmd = ["python",SCRIPTS_HOME+"detect_chimera_from_blastx_modifed.py",blastx_output_file,DIR]
		print((" ".join(cmd)))
		os.system(" ".join(cmd))
	
	assert os.path.exists(DIR+cut_file) and os.path.exists(DIR+info_file), "detect_chimera_from_blastx did not finish"
	

#Remove chimera transcripts from original file - returns non_chimera and chimera fasta files
def remove_chimeras_from_fasta(transcripts,info_file,DIR):

	if DIR == ".": DIR = os.getcwd()	
	if os.path.isabs(DIR) == False: DIR = os.path.abspath(DIR)
	if DIR[-1] != "/": DIR += "/"

	path_transcript, file_transcript = os.path.split(transcripts) #splits the path from the file name
	transcript_name = str(file_transcript)
	transcript_base_name = transcript_name.split( "." )
	
	filtered_transcripts = transcript_base_name[0]+".filtered_transcripts.fa"
	chimeric_transcripts = transcript_base_name[0]+".chimera_transcripts.fa"
	
	if os.path.exists(DIR+filtered_transcripts):
		print("Filtered transcripts file found for", transcript_base_name[0])
		
	else:
		
		transcripts_original = SeqIO.index(transcripts, "fasta")

		df = pd.read_table(info_file, header=None)
		blastx_columns = "qseqid qlen sseqid slen frames pident nident length mismatch gapopen qstart qend sstart send evalue bitscore empty".strip().split(' ')
		df.columns = blastx_columns
		chimera_names =  (df["qseqid"]).drop_duplicates()

		if len(chimera_names) == 0:
			print("No chimeras found")
		else:
			chimeras = []
			for i in chimera_names: 
			    chimeras.append(transcripts_original[i])
			count = SeqIO.write(chimeras,(DIR+chimeric_transcripts), "fasta")
			print(("Removed %i chimeras" % count))
	
			chimera_list = chimera_names.tolist()
			all_transcripts_ids = sorted(transcripts_original)
			non_chimeras_names = [x for x in all_transcripts_ids if x not in chimera_list]

			assert os.path.exists(DIR+chimeric_transcripts), "Chimera filtering did not finish"
		
		if len(non_chimeras_names) == 0:
		    print("No non-chimeras found")
		else:
			non_chimeras =[]
			for i in non_chimeras_names:
			    non_chimeras.append(transcripts_original[i])
			count = SeqIO.write(non_chimeras,(DIR+filtered_transcripts), "fasta")
			print(("Retained %i transcripts" % count))
		
			assert os.path.exists(DIR+filtered_transcripts), "Chimera filtering did not finish"
			

if __name__ == "__main__":
	if len(sys.argv) != 5:
		print ("Usage:")
		print ("python run_chimera_detection transcripts_fasta reference_proteome_fasta num_cores output_dir")
		sys.exit()

	transcripts=sys.argv[1]
	proteome_ref=sys.argv[2]
	num_cores=int(sys.argv[3])
	DIR=sys.argv[4]

	if DIR == ".": DIR = os.getcwd()	
	if os.path.isabs(DIR) == False: DIR = os.path.abspath(DIR)
	if DIR[-1] != "/": DIR += "/"

	path_ref, file_ref = os.path.split(proteome_ref) #splits the path from the file name
	ref_name = str(file_ref)
	ref_base_name = ref_name.split( "." )

	path_transcript, file_transcript = os.path.split(transcripts) #splits the path from the file name
	transcript_name = str(file_transcript)
	transcript_base_name = transcript_name.split( "." )
	blastx_output_file = transcript_base_name[0]+".blastx"

	path_blastx, file_blastx = os.path.split(DIR+blastx_output_file) #splits the path from the file name
	blastx_name = str(file_blastx)
	blastx_base_name = blastx_name.split( "." )
	info_file = (blastx_base_name[0]+".info")

	
	make_blast_db(proteome_ref,DIR)
	
	run_blastx(transcripts,((DIR+ref_base_name[0])),num_cores,DIR)
	
	run_chimera_detection((DIR+blastx_output_file),DIR)
	
	remove_chimeras_from_fasta(transcripts,(DIR+info_file),DIR)
