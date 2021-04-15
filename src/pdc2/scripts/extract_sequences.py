"""
Extract sequences from Chloroplast and Mitochondrial databases given a plant Order name (it also works for family, genus or species names and genbank accession)
"""

import sys, os, re
from Bio import SeqIO

from os.path import expanduser
home = expanduser("~")

CP_DATABASE=os.path.expanduser("~/Desktop/botany_2018/databases/chloroplast_NCBI_reference_sequences_OCT_17_2018.fasta.bgz")
MT_DATABASE=os.path.expanduser("~/Desktop/botany_2018/databases/mitochondrion_NCBI_reference_sequences_OCT_17_2018.fasta.bgz")

def extract_order_cp(order,DIR):
	assert os.path.exists(CP_DATABASE),"Cannot find the plastome database "+CP_DATABASE

	if DIR == ".": DIR = os.getcwd()	
	if os.path.isabs(DIR) == False: DIR = os.path.abspath(DIR)
	if DIR[-1] != "/": DIR += "/"
	
	extracted_order_cp_seqs = order+"_cp.fa"
	
	if os.path.exists(DIR+extracted_order_cp_seqs):
		print "Found", extracted_order_cp_seqs
	else:
		cp_db = SeqIO.index(CP_DATABASE, "fasta") #imports the bgz compressed chloroplast database file as a dictionary-like (not keept in memory)
		cp_db_keys=list(cp_db.keys()) #create a list of dictionaries keys, in this case the sequence header
		order_search=re.compile(".*"+str(order)+".*") #compile the pattern to search in the key list, in this case will be the plant Order
		order_keys = filter(order_search.match, cp_db_keys) #make a list on only the keys of the specified plant Order
		if len(order_keys) == 0: #counts key matches
			print order, "not found in chloroplast database"
		else:
			sequences = [] #opens a empty list to append the sequences of the specified order
			for i in order_keys: sequences.append(cp_db[i]) #search in indexed cp database sequences with the specified plant Order and appends to the list
			count = SeqIO.write(sequences, (DIR+order+"_cp.fa"), "fasta") #counts and writes the fasta file containing only the sequences from the specified plant Order
			print("Extracted %i chloroplast sequences for" % count),order,"in",(order+"_cp.fa") #print the number of sequences extracted and the name of the destination file
	#assert os.path.exists(DIR+extracted_order_cp_seqs), "No choloroplast file with " + order + " sequences was created"



def extract_order_mt(order,DIR):
	assert os.path.exists(MT_DATABASE),"Cannot find the mitochondrial database "+MT_DATABASE

	if DIR == ".": DIR = os.getcwd()	
	if os.path.isabs(DIR) == False: DIR = os.path.abspath(DIR)
	if DIR[-1] != "/": DIR += "/"

		
	extracted_order_mt_seqs = order+"_mt.fa"
	
	if os.path.exists(DIR+extracted_order_mt_seqs):
		print "Found", extracted_order_mt_seqs
	else:
		mt_db = SeqIO.index(MT_DATABASE, "fasta")
		mt_db_keys=list(mt_db.keys()) 
		order_search=re.compile(".*"+str(order)+".*")
		order_keys = filter(order_search.match, mt_db_keys)
		if len(order_keys) == 0: 
			print order, "not found in mitochondrial database"
		else:
			sequences = [] 
			for i in order_keys: sequences.append(mt_db[i])
			count = SeqIO.write(sequences, (DIR+order+"_mt.fa"), "fasta")
			print("Extracted %i mitochondrial sequences for" % count),order,"in",(order+"_mt.fa")
	#assert os.path.exists(DIR+extracted_order_mt_seqs), "No mitochondrial file with " + order + " sequences was created"


def extract_both_cat(order,DIR):
	
	if DIR == ".": DIR = os.getcwd()	
	if os.path.isabs(DIR) == False: DIR = os.path.abspath(DIR)
	if DIR[-1] != "/": DIR += "/"


	extract_order_cp(order,DIR)
	extract_order_mt(order,DIR)
		
	extracted_order_cp_seqs = order+"_cp.fa"
	extracted_order_mt_seqs = order+"_mt.fa"
	extracted_order_cp_mt_seqs = order+"_cp_mt.fa"

	if os.path.exists(DIR+extracted_order_cp_mt_seqs):
		print "Found", extracted_order_cp_mt_seqs
	else:	
		if os.path.exists(DIR+extracted_order_cp_seqs) and os.path.exists(DIR+extracted_order_mt_seqs):	
			extracted_sequences = [DIR+extracted_order_cp_seqs, DIR+extracted_order_mt_seqs]	
			with open((DIR+extracted_order_cp_mt_seqs), 'w') as outfile:
				for files in extracted_sequences:
					with open(files) as infile:
						outfile.write(infile.read())
			print "Concatenated chloroplast and mitochondrial file for",order,"saved as",(order+"_cp_mt.fa")		
		elif os.path.exists(DIR+extracted_order_cp_seqs) and not os.path.exists(DIR+extracted_order_mt_seqs):
			print "Only chloroplast sequences found for",order
		elif os.path.exists(DIR+extracted_order_mt_seqs) and not os.path.exists(DIR+extracted_order_cp_seqs):
			print "Only mitochondrial sequences found",order
	
	#assert os.path.exists(DIR+extracted_order_cp_mt_seqs), "No contatenated file with " + order + " sequences was created"


if __name__ == "__main__":
	if len(sys.argv) != 4:
		print "Usage:"
		print "For sequence extraction: python extract_sequences.py Order_name genome[cp, mt or both] output_dir "
	else:	
		assert sys.argv[2] == "cp" or sys.argv[2] == "mt" or sys.argv[2] == "both", \
		"genome type has to be either cp or mt or both"
		if sys.argv[2] == "cp":
			extract_order_cp(order=sys.argv[1],DIR=sys.argv[3])
		elif sys.argv[2] == "mt":
			extract_order_mt(order=sys.argv[1],DIR=sys.argv[3])
		elif sys.argv[2] == "both":
			extract_both_cat(order=sys.argv[1],DIR=sys.argv[3])
	sys.exit(0)		



