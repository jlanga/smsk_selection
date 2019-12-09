"""
Read the concatenated fasta file either with or without ends cut
write individual fasta files for each cluster
"""

import sys,os
#from Bio import SeqIO
from seq import read_fasta_file

def mcl_to_fasta(all_fasta,mcl_outfile,minimal_taxa,outdir):
	print "Reading mcl output file"
	clusterDICT = {} #key is seqID, value is clusterID
	count = 0
	if outdir[-1] != "/": outdir += "/"
	with open(mcl_outfile,"rU") as infile:
		for line in infile:
			if len(line) < 3: continue #ignore empty lines
			spls = line.strip().split('\t')
			if len(set(i.split("@")[0] for i in spls)) >= minimal_taxa:
				count += 1
				clusterID = str(count)
				for seqID in spls:
					clusterDICT[seqID] = clusterID
	print count,"clusters with at least",minimal_taxa,"taxa read"
					
	print "Reading the fasta file",all_fasta
	#handle = open(all_fasta,"rU")
	#for record in SeqIO.parse(handle,"fasta"):
	for s in read_fasta_file(all_fasta):
		#seqid,seq = str(record.id),str(record.seq)
		seqid,seq = s.name,s.seq
		try:
			clusterID = clusterDICT[seqid]
			with open(outdir+"cluster"+clusterID+".fa","a") as outfile:
				outfile.write(">"+seqid+"\n"+seq+"\n")
		except:
			pass # Those seqs that did not go in a cluster with enough taxa
			# will not be in clusterDICT
	#handle.close()
	

if __name__ =="__main__":
	if len(sys.argv) != 5:
		print "usage: write_fasta_files_from_mcl.py all_fasta mcl_outfile minimal_taxa outDIR"
		sys.exit()
	
	mcl_to_fasta(all_fasta=sys.argv[1],mcl_outfile=sys.argv[2],\
				 minimal_taxa=int(sys.argv[3]),outdir=sys.argv[4])
	
