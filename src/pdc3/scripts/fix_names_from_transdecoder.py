"""
fix EST names to taxonID@seqID before transdecoder
"""

import os,sys

#seq names look like
#>TR10000|c0_g1_i1|m.14496 TR10000|c0_g1_i1|g.14496  ORF TR10000|c0_g1_i1|g.14496 TR10000|c0_g1_i1|m.14496 type:5prime_partial len:180 (+) TR10000|c0_g1_i1:1-540(+)
#>cds.Mecr@gi_11549800_gb_BF478973_1_BF478973|m.9518 Mecr@gi_11549800_gb_BF478973_1_BF478973|g.9518  ORF Mecr@gi_11549800_gb_BF478973_1_BF478973|g.9518 Mecr@gi_11549800_gb_BF478973_1_BF478973|m.9518 type:3prime_partial len:195 (+) Mecr@gi_11549800_gb_BF478973_1_BF478973:86-667(+)

if __name__ =="__main__":
	if len(sys.argv) != 3:
		print("usage: python fix_names_from_transdecoder.py inDIR outDIR")
		sys.exit()
	
	inDIR = sys.argv[1]+"/"
	outDIR = sys.argv[2]+"/"
	for i in os.listdir(inDIR):
		taxonID = i.split(".")[0]
		if i.endswith("transdecoder.pep"):
			outname = taxonID+".pep.fa"
		elif i.endswith("transdecoder.cds"):
			outname = taxonID+".cds.fa"
		elif i.endswith(".Trinity.fasta"):
			os.system("gzip "+inDIR+i)
			continue
		else: continue
		print(i)
		
		infile = open(inDIR+i,"rU")
		outfile = open(outDIR+outname,"w")
		for line in infile:
			if line[0] == ">":
				newid = line.split(" ")[0]
				newid = newid.split(".")[-1]
				outfile.write(">"+taxonID+"@"+newid+"\n")
			else: outfile.write(line)
		outfile.close
		infile.close()
				
		os.system("gzip "+inDIR+i)
	
