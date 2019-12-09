"""
The blastx output files must have a customized tubular output:
0-qseqid	1-qlen		2-seqid		3-slen		4-frame		5-pident
6-nident	7-length	8-mismatch	9-gapopen	10-qstart	11-qend
12-start	13-send		14-evalue	15-bitscore

Step 1: only look at HSPs that are long and sufficiently similar to the query
Step 2: check a block of HSPs from the a query-hit pair (hit_block) for self chimeras
Step 3: check  a block of HSPs from the same query (query_block) for multi-gene chimeras
		if no self chimera was detected. Out put two cut ranges with a and b after the names

output a ".cut" file for position to cut, and a ".info" file for validating cuts

Modify the PIDENT_CUTOFF and LENGTH_CUTOFF as needed for each data set
The difference between "trans-multi" vs. "trans-self" is only meaningful when
the blast database is from a single species
"""

import sys,os

PIDENT_CUTOFF = 30 #only look at HSPs >= this percentage similarity
LENGTH_CUTOFF = 100	#only look at HSPs >= this length

#calculate length of query coverage
def qcov(hsp):
	return abs(hsp[11]-hsp[10]) + 1
	
#given two hsps, return True if 
#overlap les than 20% of the shorter and overlap less than 60 bp
def separated(hsp1,hsp2):
	length1 = qcov(hsp1)
	length2 = qcov(hsp2)
	start = min(hsp1[10],hsp1[11],hsp2[10],hsp2[11])
	end = max(hsp1[10],hsp1[11],hsp2[10],hsp2[11])
	overlap = length1+length2 - (end-start) + 1
	#value of overlap can < 0 but only the upper limit maters
	if overlap < min(60, 0.2*min(length1,length2)):
		return True
	else: return False


#expand query range given two hsps of the same query-hit pair
#both hsps are the same direction
def expand_range(hsp1,hsp2):
	if hsp1 == []: return hsp2
	if hsp2 == []: return hsp1
	start1,end1,start2,end2 = hsp1[10],hsp1[11],hsp2[10],hsp2[11]
	if start1 < end1 and start2 < end2:#both forward
		start,end = min(start1,start2),max(end1,end2)
	elif start1 > end1 and start2 > end2:#both reverse
		start,end = max(start1,start2),min(end1,end2)
	#no change if new hsp is of opposite direction
	hsp1[10],hsp1[11] = start,end
	return hsp1
	

#detect chimera from a block of hits
#block can be a query-hit block (only check for self-chimera)
#or block can be a query block (check both self and multi-gene chimera)
#return True if no chimera is detected, False if chimera is detected
#also write to out file the combined best HSP
def check_block(block,multigene):
	#only one hsp -> not chimera
	if len(block) == 1: return True	
	#summarize all pos and neg HSPs
	pos,neg = [],[]
	for hsp in block:
		if hsp[4][0] == "-":
			neg = expand_range(neg,hsp)
		else:
			pos = expand_range(pos,hsp)
	#compare pos_hit and neg_hit
	if (pos == [] and neg != []) or (neg == [] and pos != []):
		return True #only has hits of one direction
	elif separated(pos,neg):#has both direction and separate -> chimera!
		#write range to cut
		if multigene: #output both hsps
			start1,end1 = min(pos[10],pos[11]),max(pos[10],pos[11])
			start2,end2 = min(neg[10],neg[11]),max(neg[10],neg[11])
			outfile1.write(pos[0]+" "+str(int(start1))+" "+str(int(end1))+" trans-multi\n")
			outfile1.write(neg[0]+" "+str(int(start2))+" "+str(int(end2))+" trans-multi\n")
		else:
			if qcov(pos) > qcov(neg):
				outhsp = pos #outhsp is the better covered of the two
			else: outhsp = neg
			start,end = min(outhsp[10],outhsp[11]),max(outhsp[10],outhsp[11]) #range to cut
			outfile1.write(outhsp[0]+" "+str(int(start))+" "+str(int(end))+" trans-self\n")
		#write the blastx block to a .info file for visual checking
		for i in pos:
			outfile2.write(str(i)+"\t")
		outfile2.write("\n")
		for i in neg:
			outfile2.write(str(i)+"\t")
		outfile2.write("\n")
		return False
	else:
		return True #has both direction but not separate


if __name__ == "__main__":
	if len(sys.argv) != 3:
		print "usage: python detect_chimera_from_blastx.py blastx_output output_dir"
		sys.exit()

	blastx_output = sys.argv[1]
	DIR = sys.argv[2]
	
	if DIR == ".": DIR = os.getcwd()	
	if os.path.isabs(DIR) == False: DIR = os.path.abspath(DIR)
	if DIR[-1] != "/": DIR += "/"

	path_blastx, file_blastx = os.path.split(blastx_output) #splits the path from the file name
	blastx_name = str(file_blastx)
	blastx_base_name = blastx_name.split( "." )

	infile = open(blastx_output,"rU")
	outfile1 = open(DIR+blastx_base_name[0]+".cut","w")
	outfile2 = open(DIR+blastx_base_name[0]+".info","w")
	last_query = ""
		
	
	for line in infile:
		if len(line) < 3: continue #ignore empty lines
		hsp = line.strip().split("\t")
		for i in [5,10,11]: hsp[i] = float(hsp[i])
		if hsp[5] < PIDENT_CUTOFF or qcov(hsp) < LENGTH_CUTOFF:
			continue #ignore low similarity or short HSPs
		query,hit = hsp[0],hsp[2]
		
		if last_query == "": #at the beginning of a file
			hit_block = [hsp] #store all HSPs of the same query and same hit
			query_block = [hsp] #store all HSPs from the same query
			good_seq = True #False if chimera is detected
		
		elif query == last_query: #continue with the same query
			query_block.append(hsp)
			if good_seq: #only check if no chimera has been detected
				if hit == last_hit:
					hit_block.append(hsp)
				else: #send off the hit_block
					good_seq = check_block(hit_block,False)
					hit_block = [hsp]

		else: #starting a new query
			if good_seq: #haven't found self chimera
				good_seq = check_block(hit_block,False) #check the last hit block
			if good_seq: #haven't found self chimera
				good_seq = check_block(query_block,True) #look for multi-chimera
			query_block,hit_block = [hsp],[hsp]
			good_seq = True
			
		#keep track of the last line processed	
		last_query,last_hit = query,hit
	
	if good_seq: #haven't found self chimera
		good_seq = check_block(hit_block,False) #check the last hit block
	if good_seq: #haven't found self chimera
		good_seq = check_block(query_block,True) #check the last query block
	
	infile.close()
	outfile1.close()
	outfile2.close()