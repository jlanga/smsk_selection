"""
copy both the .phy and .model file into the same working dir
cd into the working dir

Make exclude files for each jackknife replicate and run it 
Also print out a summary file with each line being the gene number that 
was excluded at each jackknife replicate
"""

import os,sys,random

REPEATS = 200 #generate 200 jackknife exclude files
RAXML_CMD = "raxml"

if __name__ == "__main__":
	if len(sys.argv) != 4:
		print "usage: python jackknife_by_percent_genes.py num_core jacknife_proportion(between 0 and 1) dna/aa"
		sys.exit()
	
	num_core = sys.argv[1]
	JK_PROP = float(sys.argv[2])
	if JK_PROP <= 0 or JK_PROP > 1:
		print "jacknife proportion has to be between 0 and 1"
		sys.exit()
		
	if sys.argv[3] == "dna":
		MODEL = "GTRCAT"
	elif sys.argv[3] == "aa":
		MODEL = "PROTCATWAG"	
	else:
		print "Input data type: DNA or aa"
		sys.exit()
		
	for i in os.listdir("./"):
		if i[-6:] == ".model": model_file = i
		elif i[-4:] == ".phy": aln = i
	
	#parse the .model file and get the partitions
	infile = open(model_file,"r")
	geneDICT = {} #key is gene id, value is (start,end)
	for line in infile:
		if len(line) < 3: continue
		#line looks like DNA,cluster476-2-1_1.subtree1.p1.aln-cln_gene1=1-1216
		spls = ((line.strip()).split("gene")[1]).split("=")
		geneid = spls[0]
		rg = spls[1].split("-") #range of position in supermatrx
		geneDICT[geneid] = int(rg[0]),int(rg[1])
		
	#randomly generate exclude files and run jackknife
	gen_num = len(geneDICT)
	print gen_num,"genes in total"
	jk_num = int((gen_num * (1-JK_PROP))+0.5)
	print jk_num,"genes excluded in each round of jackknife analysis"
	outfile1 = open("jacknife_summary","w")
	outfile1.write("total number of genes "+str(gen_num)+"\n"+str(jk_num)+"excluded for each replicate\n")
	JK_PROP = str(int(JK_PROP*100))
	for i in range(REPEATS):
		exclude = [] #geneIDs to exclude
		while True:
			r = random.randint(1,gen_num)
			if r not in exclude:
				exclude.append(r)
			if len(exclude) == jk_num:
				break
		exclude.sort()
		exc_file = "jk"+JK_PROP+"_rep"+str(i+1)
		with open(exc_file,"w") as outfile:
			for j in exclude:
				outfile.write(str(geneDICT[str(j)][0])+"-"+str(geneDICT[str(j)][1])+" ")
				outfile1.write(str(j)+" ")
		outfile1.write("\n")
		#run raxml for this jackknife replicate
		cmd = RAXML_CMD+" -F -T "+num_core+" -E "+exc_file+" -m "+MODEL+" -p 6666 -q "+model_file
		cmd += " -s "+aln+" -n "+aln.replace("phy",exc_file)
		print cmd
		os.system(cmd) #write new matrix
		cmd = RAXML_CMD+" -F -T "+num_core+" -m "+MODEL+" -p 6666 -q "+model_file+"."+exc_file
		cmd += " -s "+aln+"."+exc_file+" -n "+aln.replace("phy",exc_file)
		print cmd
		os.system(cmd) #run raxml on the subsampled alignment
	outfile1.close()
		
