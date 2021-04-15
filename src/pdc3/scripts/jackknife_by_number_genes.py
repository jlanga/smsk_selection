"""make sure that raxml is in the path"""

import os,sys,random
import subprocess

def jk(indir,num_core,resample_num,seqtype,replicates=200):
	"""
	Write an exclude files for each jackknife replicate
	Also print out a summary file with each line being the gene number that 
	was excluded at each jackknife replicate
	"""
	if indir[-1] != "/": indir += "/"
	outdir = indir+str(resample_num)+"genes/"
	try: os.stat(outdir)
	except: os.mkdir(outdir)
	count_model, count_phy = 0,0
	for i in os.listdir(indir):
		if i.endswith(".model"):
			model_file = i if indir == "./" else indir+i
			count_model += 1
		elif i.endswith(".phy"): 
			aln = i if indir == "./" else indir+i
			count_phy += 1
	assert count_model == 1 and count_phy == 1,\
		"Check to make sure that only 1 .model and 1 .phy file in indir"
	assert seqtype == "dna" or seqtype == "aa", \
		"seqtype: dna or aa"
	model = "GTRCAT" if seqtype == "dna" else "PROTCATWAG"

	print("Found sequence alignment",aln)
	print("Found partition file",model_file)
	print("Resample",resample_num,"each time for",replicates,"replicates")
	
	# parse the .model file and get the partitions
	geneDICT = {} #key is gene id, value is range of columns in a string
	with open(model_file,"r") as infile:
		for line in infile:
			if len(line) < 3: continue
			# line looks like DNA,cluster2772rr.1to1ortho.aln-cln_gene22=23371-24283
			spls = ((line.strip()).split("gene")[1]).split("=")
			# spls look like ["22","23371-24283"]
			geneDICT[spls[0]] = spls[1]
		
	#randomly generate exclude files and run jackknife
	gen_num = len(geneDICT)
	print(gen_num,"genes in total")
	jk_num_to_exclude = gen_num - int(resample_num)
	print(jk_num_to_exclude,"genes excluded in each jackknife replicate")
	outfile1 = open(outdir+"jacknife_summary","w")
	outfile1.write("total number of genes "+str(gen_num)+"\n")
	outfile1.write(str(jk_num_to_exclude)+" excluded for each replicate\n")
	for i in range(replicates):
		exclude = [] #geneIDs to exclude
		while len(exclude) < jk_num_to_exclude:
			r = random.randint(1,gen_num) # [1,gen_num]
			if r not in exclude:
				exclude.append(r)
		exclude.sort()
		repid = str(resample_num)+"genes_rep"+str(i+1)
		with open(repid,"w") as outfile:
			for j in exclude:
				geneID = str(j)
				outfile.write(geneDICT[geneID]+"\n")
				outfile1.write(geneID+" ")
		outfile1.write("\t")
		
		# write matrix for the current jackknife replicate
		cmd = ["raxml",\
			   "-E", repid,\
			   "-T", str(num_core),\
			   "-m", model,\
			   "-p", "6666",\
			   "-q", model_file,\
			   "-s", aln,\
			   "-n", repid]
		print(" ".join(cmd))
		p = subprocess.Popen(cmd,stdout=subprocess.PIPE)
		out = p.communicate()
		assert p.returncode == 0,"Error raxml"+out[0]
		
		# run raxml
		cmd = ["raxml","-F",\
			   "-T", str(num_core),\
			   "-m", model,\
			   "-p", "6666",\
			   "-q", model_file+"."+repid,\
			   "-s", aln+"."+repid,\
			   "-n", repid]
		print(" ".join(cmd))
		p = subprocess.Popen(cmd,stdout=subprocess.PIPE)
		out = p.communicate()
		assert p.returncode == 0,"Error raxml"+out[0]

		try: # remove intermediate files
			os.rename("RAxML_result."+repid,outdir+"RAxML_result."+repid)
			os.rename(repid,outdir+repid)
			os.remove("RAxML_info."+repid)
			os.remove("RAxML_log."+repid)
			os.remove("RAxML_parsimonyTree."+repid)
			os.remove(model_file+"."+repid)
			os.remove(aln+"."+repid)
			os.remove(aln+"."+repid+".reduced")
		except: pass # no need to worry about extra intermediate files
		
	os.system("cat "+outdir+"RAxML_result.* >"+str(resample_num)+"genes_trees")
	outfile1.close()
			

if __name__ == "__main__":
	if len(sys.argv) != 5:
		print("usage: python jackknife_by_number_genes.py indir num_core num_genes_to_resample dna/aa")
		sys.exit()
		
	indir,num_core,resample_num,seqtype = sys.argv[1:]
	jk(indir,num_core,resample_num,seqtype)
