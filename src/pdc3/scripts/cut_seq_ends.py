"""
input:
- a fasta file with all sequences used for all-by-all blast
- a file with all the unfiltered results from all-by-all blastn or blastp
Currently assume that query and hit are the same direction

Ignore hits from the same taxa
Check for ends that doesn't have any hits in any other taxa

output:
- .cutinfo file with start and end of seq to keep, and size of seq
- .cut file with all the sequences after cutting
"""

import os,sys
from Bio import SeqIO

MIN_SEQ_LEN = 40

#if taxon id pattern changes, change it here
def get_taxonid(seqid):
	return seqid.split("@")[0]

def cut_seq_ends(fasta,blast_output,logfile="log"):
	print("Reading raw blast output")	
	cutDICT = {} #key is seqid, value is a list [start,end,length]
	with open(blast_output,"r") as infile:
		for line in infile:
			if len(line) < 3:
				continue #skip empty lines
			spls = line.strip().split("\t")
			query,hit = spls[0],spls[2]
			if get_taxonid(query) == get_taxonid(hit):
				continue #ignore hits from the same taxa
			qlen,qstart,qend = int(spls[1]),int(spls[10]),int(spls[11])
			slen,sstart,send = int(spls[3]),int(spls[12]),int(spls[13])

			#get the widest range
			if query not in cutDICT:
				cutDICT[query] = [10000000,1,qlen] #[start,end,qlen]
			if hit not in cutDICT:
				cutDICT[hit] = [10000000,1,slen] #[start,end,slen]
			cutDICT[query][0] = min(cutDICT[query][0],qstart,qend) #compare starts
			cutDICT[query][1] = max(cutDICT[query][1],qstart,qend) #compare ends
			cutDICT[hit][0] = min(cutDICT[hit][0],sstart,send) #compare starts
			cutDICT[hit][1] = max(cutDICT[hit][1],sstart,send) #compare ends

	#output seqid, start and end for cutting, and seq length
	with open(blast_output+".cutinfo","w") as outfile:
		for seqid in cutDICT:
			start,end,length = cutDICT[seqid] #[start,end,length]
			outfile.write(seqid+"\t"+str(start)+"\t"+str(end)+"\t"+str(length)+"\n")
	print("Output written to",sys.argv[1]+".cutinfo")

	print("Cutting")
	outfile = open(fasta+".cut","w")
	outdict = {} # key is taxonid, value is [before,after,half_left]
	with open(fasta,"r") as handle:
		for record in SeqIO.parse(handle,"fasta"):
			seqid,seq = str(record.id),str(record.seq)
			taxonid = get_taxonid(seqid)
			if taxonid not in outdict:
				outdict[taxonid]= [0,0,0]
			outdict[taxonid][0] += 1
			if seqid in cutDICT:
				start,end,length = cutDICT[seqid]
				seq_cut = seq[start-1:end]
				if len(seq_cut) >= (len(seq)/2):
					outdict[taxonid][1] += 1 # at least half survived cutting
				#print seqid, start, end, length,end-start+1,MIN_SEQ_LEN
				if len(seq_cut) >= MIN_SEQ_LEN:
					outfile.write(">"+seqid+"\n"+seq_cut+"\n")
					outdict[taxonid][2] += 1
			else: pass # remove seqs with no interspecific hits
	outfile.close()
	summary = "taxonid\tbefore\thalf_left\tafter\tperc_survive\n"
	for i in outdict:
		out = outdict[i]
		summary+= i+"\t"+str(out[0])+"\t"+str(out[1])+"\t"+str(out[2])+"\t"
		summary+= str(out[2]/float(out[0])*100)+"%\n"
	print(summary)
	with open(logfile,"a") as f: f.write(summary)
	return fasta+".cut"

if __name__ =="__main__":
	if len(sys.argv) != 3:
		print("usage: cut_seq_ends.py fasta blast_output")
		sys.exit()
	
	fasta,blast_output = sys.argv[1:]
	cut_seq_ends(fasta,blast_output)
	

