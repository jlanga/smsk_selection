"""
input: concatenat all the raw all-by-all blast results into one file

blast output file columns (separated by tab):
0-qseqid 1-qlen 2-sseqid 3-slen 4-frames 5-pident 6-nident 7-length 
8-mismatch 9-gapopen 10-qstart 11-qend 12-sstart 13-send 14-evalue 15-bitscore

Ignore self-hits, filter by hit fraction, and convert evalues to -log(Evalue)
If want to exclude certain taxa, put the taxon IDs to exclude in EXCLUDE

if -logEvalue < 0.01, output 0.01
just arbitorily pick a really small vallue to avoid netagive values

Check for possible contamination between two data sets and output hits that
are 100% identical and overlap for > 90% of the shorter one
"""

import os,sys
import math

#EXCLUDE = ["Pinu","Pice","Zmay","Mesc","Tcac"] #taxon IDs to ignore
EXCLUDE = []
IGNORE_INTRASPECIFIC_HITS = False

def get_taxon_name(name):
	return name.split("@")[0]

def get_minus_log_evalue(raw_value):
	try:
		result = -math.log10(float(raw_value))
	except:
		result = 180.0 #largest -log10(evalue) seen
	 # mcl cannot take negative values
	 # have to use a value close to zero for all the values that are lower
	if result < 0.01:
		result = 0.01
	return result

def blast_to_mcl(blast_output,hit_frac_cutoff,exclude=EXCLUDE,\
				 ignore_interspecific_hits=IGNORE_INTRASPECIFIC_HITS):

	#both percent query and hit coverage have to >= this
	
	print("Reading raw blast output")
	print("Taxa to exclude:", EXCLUDE)
	print("Ignore intraspecific hit?",IGNORE_INTRASPECIFIC_HITS)
	infile = open(blast_output,"r")
	outname = blast_output+".hit-frac"+str(hit_frac_cutoff)+".minusLogEvalue"
	if len(EXCLUDE) > 0:
		outname += ".exclude"
	if IGNORE_INTRASPECIFIC_HITS:
		outname += ".interspecific"
	print("Writing output to",outname)
	outfile = open(outname,"w")
	count = 0
	last_query,last_hit = "",""
	print("Warning: highly similar sequences between data sets.")
	print("written to file",blast_output+".ident")
	outfile1 = open(blast_output+".ident","w")
	for line in infile:
		if len(line) < 3: continue #skip empty lines
		spls = line.strip().split("\t")
		query,hit = spls[0],spls[2]
		if query == hit:
			continue #skip self hits
		query_taxon, hit_taxon = get_taxon_name(query),get_taxon_name(hit)
		if query_taxon in EXCLUDE or hit_taxon in EXCLUDE:
			continue
		if IGNORE_INTRASPECIFIC_HITS and query_taxon == hit_taxon:
			continue
		pident = float(spls[5])
		qlen,qstart,qend = float(spls[1]),float(spls[10]),float(spls[11])
		slen,sstart,send = float(spls[3]),float(spls[12]),float(spls[13])
		
		# Check for contamination
		if pident == 100.0 and query_taxon != hit_taxon:
			if abs(qend-qstart)/qlen >= 0.9 or abs(send-sstart)/slen >= 0.9:
				if qlen >= 300 and slen >= 300: # ignore short seqs
					print(line, end=' ')
					outfile1.write(line)
		
		minusLogEvalue = get_minus_log_evalue(spls[14])
		if last_query == "":
			#The highest -log(evalue) for each query-hit pair that is not a self hit
			#is the first one encoutered in the blast output file
			max_minusLogEvalue = minusLogEvalue
		else: #not at the very beginning of the blastn output
			if query == last_query and hit == last_hit:
				#expand the query range and hit ranges
				#all the query and hit are in the same direction for aa and CDS
				qstart = min(qstart,last_qstart)
				qend   = max(qend,last_qend)
				sstart = min(sstart,last_sstart)
				send   = max(send,last_send)
			else:#summarize last query-hit pair
				perc_qrange = (last_qend - last_qstart + 1) / last_qlen
				perc_srange = (last_send - last_sstart + 1) / last_slen
				if perc_qrange >= float(hit_frac_cutoff) and perc_srange >= float(hit_frac_cutoff):
					#output info for the previous query-hit pair
					outfile.write(last_query+"\t"+last_hit+"\t"+str(max_minusLogEvalue)+"\n")
					count += 1
					if count % 1000000 == 0:
						print(str(count/1000000),"million hits written to the output")
				max_minusLogEvalue = minusLogEvalue #reset max_minusLogEvalue for a new query-hit pair
		last_query,last_qlen,last_qstart,last_qend = query,qlen,qstart,qend
		last_hit,last_slen,last_sstart,last_send = hit,slen,sstart,send
	
	# process the last
	perc_qrange = (last_qend - last_qstart + 1) / last_qlen
	perc_srange = (last_send - last_sstart + 1) / last_slen
	if perc_qrange >= float(hit_frac_cutoff) and perc_srange >= float(hit_frac_cutoff):
		outfile.write(last_query+"\t"+last_hit+"\t"+str(max_minusLogEvalue)+"\n")
	infile.close()
	outfile.close()
	outfile1.close()
	print("Output written to",outname)
	print("Highly similar interspecific hits written to",blast_output+".ident")
	print("Check for possible contamination")
	return outname
			

if __name__ == "__main__":
	if len(sys.argv) != 3:
		print("usage: python blast_to_mcl.py rawblast hit_fraction")
		sys.exit()
	
	blast_to_mcl(blast_output=sys.argv[1],hit_frac_cutoff=float(sys.argv[2]))
