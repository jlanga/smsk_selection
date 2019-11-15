"""
Check data set names for:
- Duplicated names
- Special characters other than digits, letters and _
- all names look like taxonID@seqID, and file names are the taxonID
"""

import os,sys
from re import match

if __name__ =="__main__":
	if len(sys.argv) != 3:
		print("usage: python check_names.py DIR file_ending")
		sys.exit()
	
	DIR = sys.argv[1]+"/"
	file_ending = sys.argv[2]
	for i in os.listdir(DIR):
		if not i.endswith(file_ending): continue
		taxonID = i.split(".")[0]
		print(i,taxonID)
		infile = open(DIR+i,"rU")
		seqIDs = []
		for line in infile:
			if line[0] != ">": continue #only look at seq names
			taxon,seqID = (line.strip()[1:]).split("@")
			assert taxon == taxonID,taxon+" taxon id do not match file name "+taxonID
			assert match("^[a-zA-Z0-9_]+$",seqID),seqID+" seq id contains special characters other than letter, digit and '_'"
			seqIDs.append(seqID)
		repeats = {j:seqIDs.count(j) for j in seqIDs}
		print(len(seqIDs),"sequences checked")
		for j in repeats:
			rep = repeats[j]
			if rep > 1:
				print(j,"repeats",rep,"times")
		infile.close()
			
		
