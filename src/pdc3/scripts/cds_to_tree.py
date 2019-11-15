import argparse
import subprocess
import shutil
import os,sys
import tree_utils
from cut_seq_ends import cut_seq_ends
from blast_to_mcl import blast_to_mcl
from write_fasta_files_from_mcl import mcl_to_fasta
import fasta_to_tree
import trim_tips
import mask_tips_by_taxonID_transcripts
import cut_long_internal_branches
import write_fasta_files_from_trees
import phyutility_wrapper
from raxml_bs_wrapper import raxml_bs
from filter_1to1_orthologs import get_121
from prune_paralogs_RT import RT
from write_alignments_from_orthologs import ortho_to_aln
from concatenate_matrices import concatenate
from jackknife_by_number_genes import jk

CDHITEST_SIMILARITY = 0.99 # cd-hit-est global similarity cutoff

def get_min_taxa(num_taxa,max_missing=1):
	""" often lots of clusters with full taxon occupancy when num_taxa is small
	don't want to print out too many clusters for a test run"""
	return num_taxa - int(max_missing) if num_taxa >= 10 else num_taxa

def check_dependencies():
	"""Make sure that all the dependencies can be located"""
	for i in ['mcl','makeblastdb','blastn','cd-hit-est','raxml',\
			   'run_pasta.py']:
		print("Checking",i,"...", end=' ')
		try:
			p = subprocess.Popen([i],stderr=subprocess.PIPE,stdout=subprocess.PIPE)
			p.communicate()
		except OSError:
			print("Could not find",i,"in the path")
			print("Make sure that it is properly installed")
			sys.exit(1)
		print("OK")
	for i in ["mafft","phyutility"]:
		print("Checking",i,"...", end=' ')
		try: os.system('which '+i)
		except: 
			print("Could not find "+i+" in the path")
			print("Make sure that it is properly installed")
			sys.exit(1)
		print("OK")

			
def gether_fasta_files(path,file_ending,logfile):
	"""find all the input files"""
	if path[:-1] != "/": path += "/"
	files,taxonIDs = [],[]
	for i in os.listdir(path):
		if i.endswith(file_ending):
			files.append(i)
			print("Reading input cds file",i)
			taxonIDs.append(i.split(".")[0])
	assert len(files) > 0, "No file ends with "+file_ending+" found in "+path
	assert len(taxonIDs) == len(set(taxonIDs)),\
		"Input files should looke like: taxonID."+file_ending
	if not os.path.exists(logfile): # only need to log input files once
		with open(logfile,"a") as f:
			for i in files: f.write(i+"\n")
	return files
	
def reduce_redundancy(indir,fasta_files,num_cores,outdir,logfile,\
		similarity=str(CDHITEST_SIMILARITY)):
	"""run cd-hit-est"""
	if indir[-1] != "/": indir += "/"
	if outdir[-1] != "/": outdir += "/"
	num_cores = str(num_cores)
	outdir1 = outdir+"1_cd-hit-est/"
	if os.path.exists(outdir+"cd-hit-est_ok"): return outdir1
	try: os.stat(outdir1)
	except: os.mkdir(outdir1)
	for i in fasta_files:
		if os.path.exists(outdir1+i+".cdhitest"): continue
		cmd = ["cd-hit-est",\
			   "-i",indir+i,\
			   "-o",outdir1+i+".cdhitest",\
			   "-c",similarity,\
			   "-n","10",
			   "-r","0",\
			   "-T",num_cores]
		print(" ".join(cmd))
		p = subprocess.Popen(cmd,stdout=subprocess.PIPE)
		out = p.communicate()
		assert p.returncode == 0,out[0]+"cd-hit-est failed"
		with open(logfile,"a") as f: f.write(" ".join(cmd)+"\n")
	os.system("touch "+outdir+"cd-hit-est_ok")
	return outdir1
	
def clustering(outdir1,fasta_files,num_cores,outdir,min_taxa,logfile):
	"""Make initial clusters by blast and mcl"""
	if outdir[-1] != "/": outdir += "/"
	outdir2 = outdir+"2_clustering/"
	try: os.stat(outdir2)
	except: os.mkdir(outdir2)
	num_cores = str(num_cores)
	
	# concatenate fasta files
	if not os.path.exists(outdir+"concatenate_all_fasta_ok"):
		cmd = ["cat"]+[outdir1+i+".cdhitest" for i in fasta_files]
		print(" ".join(cmd))
		out = open(outdir2+"all.fa", 'w')
		p = subprocess.Popen(cmd,stdout=out)
		out.close()
		p.communicate()
		assert p.returncode == 0,"Error concatenating all the cd-hit-est output"
		with open(logfile,"a") as f: f.write(" ".join(cmd)+"\n")
		os.system("touch "+outdir+"concatenate_all_fasta_ok")
		
	# make a blast database
	if not os.path.exists(outdir+"blast_database_ok"):
		cmd = ["makeblastdb",\
			   "-in",outdir2+"all.fa",\
			   "-parse_seqids",\
			   "-dbtype","nucl",\
			   "-out",outdir2+"all.fa"]
		print(" ".join(cmd))
		p = subprocess.Popen(cmd)
		p.communicate()
		assert p.returncode == 0,"Error making blast database"
		with open(logfile,"a") as f: f.write(" ".join(cmd)+"\n")
		os.system("touch "+outdir+"blast_database_ok")
	
	# blastn
	if not os.path.exists(outdir+"blastn_ok"):
		fmtstring = "'6 qseqid qlen sseqid slen frames pident nident length "
		fmtstring += "mismatch gapopen qstart qend sstart send evalue bitscore'"
		cmd =["blastn",\
			  "-db",outdir2+"all.fa",\
			  "-query",outdir2+"all.fa",\
			  "-evalue","10",\
			  "-num_threads",num_cores,\
			  "-max_target_seqs","1000",\
			  "-out",outdir2+"all.rawblast",\
			  "-outfmt",fmtstring]
		print(" ".join(cmd))
		os.system(" ".join(cmd))
		os.system("touch "+outdir+"blastn_ok")
		with open(logfile,"a") as f: f.write(" ".join(cmd)+"\n")
		os.system("rm "+outdir2+"all.fa.n*") # remove the blast database
	
	# cut seq ends that do not have any interspecific blastn hits
	if not os.path.exists(outdir+"cut_seq_ends_ok"):
		cut_seq_ends(fasta=outdir2+"all.fa",\
					 blast_output=outdir2+"all.rawblast",\
					 logfile=logfile)
		os.system("touch "+outdir+"cut_seq_ends_ok")

	# Filter raw blast output by hit fraction and mcl
	mcl_out = outdir2+"hit-frac0.4_I1.4_e5"
	if not os.path.exists(mcl_out):
		# hit-fraction of 0.4 generally works well
		mcl_infile = blast_to_mcl(blast_output=outdir2+"all.rawblast",\
								  hit_frac_cutoff=0.4)
		# inflation vallue of 1.4 generally works well for cds
		cmd = ["mcl",mcl_infile,"--abc","-te",num_cores,"-tf","'gq(5)'",\
			   "-I","1.4","-o",mcl_out]
		print(" ".join(cmd))
		os.system(" ".join(cmd))
		assert os.path.exists(outdir2+"hit-frac0.4_I1.4_e5"),"Error mcl"
		with open(logfile,"a") as f: f.write(" ".join(cmd)+"\n")
	
	# Write fasta files for each cluster from mcl output.
	outdir3 = outdir+"3_clusters/"
	try: os.stat(outdir3)
	except: os.mkdir(outdir3)
	if not os.path.exists(outdir+"mcl_to_fasta_ok"):
		mcl_to_fasta(all_fasta=outdir2+"all.fa.cut",\
					 mcl_outfile=outdir2+"hit-frac0.4_I1.4_e5",\
					 minimal_taxa=min_taxa,\
					 outdir=outdir3)
		os.system("touch "+outdir+"mcl_to_fasta_ok")
	return outdir2,outdir3

def refine(curdir,nextdir,prefix,num_cores,\
		relative_cut,absolute_cut,branch_len_cutoff,\
		min_taxa,mask_para,logfile,test=False):
	fasta_to_tree.main(curdir,num_cores,"dna",bs="n",test=test)
	trim_tips.main(DIR=curdir,\
			tree_file_ending=".tre",\
			relative_cut=relative_cut,\
			absolute_cut=absolute_cut)
	mask_tips_by_taxonID_transcripts.main(treDIR=curdir,\
			clnDIR=curdir,\
			para=mask_para,\
			intree_file_ending=".tt")
	try: os.stat(nextdir)
	except: os.mkdir(nextdir)
	cut_long_internal_branches.main(inDIR=curdir,\
		file_ending=".mm",\
		branch_len_cutoff=branch_len_cutoff,\
		min_taxa=min_taxa,\
		outDIR=nextdir)
	log = "refine in "+curdir+" with\nrelative tip cutoff of "+str(relative_cut)
	log += "\nabsolute tip cutoff of "+str(absolute_cut)
	log += "\ninternal branch cutoff of "+str(branch_len_cutoff)+"\n"
	with open(logfile,"a") as outfile: outfile.write(log)


def get_args():
	"""Get arguments """
	parser = argparse.ArgumentParser(
			description="""Example:\n
			python cds_to_homologs.py --indir=data/ --outdir=./ --threads=4 --reltip=0.2 --abstip=0.4 --deep=0.3 --ortho=RT --inout=../in_out""")
	parser.add_argument("--indir",
			help="The directory of all the taxonID.cds.fa files")
	parser.add_argument("--outdir",
			help="Output directory")
	parser.add_argument("--threads",type=int,default=4,
			help="Number of threads. Optimal 4 to 10. Default: 4")
	parser.add_argument("--reltip",type=float,default=0.2,
			help="Remove tips longer than this and > 10x of its sister.Default: 0.2")
	parser.add_argument("--abstip",type=float,default=0.4,
			help="""Remove tips longer than this. Default: 0.4""")
	parser.add_argument("--deep",type=float,default=0.4,
			help="Cut internal branches longer than this to remove deep paralogs. Default: 0.4")
	parser.add_argument("--ortho",help="121/RT")
	parser.add_argument("--inout",
			help="A text file that each line looks like IN	taxonID or OUT	taxonID. Needed if using RT for orthology inference")
	parser.add_argument("--max_mis_taxa",type=int,default=1,
			help="""Maximum number of taxa missing in homolog and orthologs. Default: 1""")
	parser.add_argument("--test", default="n",
			help="""If test is 'y', only analyze clusters end with zero""")
	return parser.parse_args()


def main():
	args = get_args()
	print(args)
	indir,outdir,num_cores = args.indir, args.outdir, args.threads
	relative_cut,absolute_cut,branch_len_cutoff=args.reltip,args.abstip,args.deep
	assert args.ortho == "RT" or args.ortho == "121", "--ortho has to be 121 or RT"
	assert os.path.exists(args.inout),"cannot find the file specified by --inout"
	assert args.test == "y" or args.test == "n", "test has to be either y or n"
	test = True if args.test == "y" else False
	if outdir[-1] != "/": outdir += "/"
	if indir[-1] != "/": indir += "/"
	
	logfile = outdir+"homology_inference.log"
	check_dependencies()
	
	# get initial fasta
	fasta_files = gether_fasta_files(path=indir,file_ending=".cds.fa",logfile=logfile)
	print(len(fasta_files),"data sets read")
	taxa = set([tree_utils.get_name(i.split(".")[0]) for i in fasta_files])
	print(len(taxa),"taxa found:")
	print(taxa)
	min_taxa = get_min_taxa(len(taxa),args.max_mis_taxa)
	print("Minimal number of taxa: ",min_taxa)
	outdir1 = reduce_redundancy(indir,fasta_files,num_cores,outdir,logfile)
	outdir2,outdir3 = clustering(outdir1,fasta_files,num_cores,outdir,min_taxa,logfile)
	
	# tree inference and clearning round 1.
	if not os.path.exists(outdir+"3_clusters_ok"):
		refine(curdir=outdir3,nextdir=outdir+"4_refine/",\
				prefix="3_clusters",\
				num_cores=num_cores,\
				relative_cut=relative_cut,\
				absolute_cut=absolute_cut,\
				branch_len_cutoff=branch_len_cutoff,\
				min_taxa=min_taxa,\
				mask_para="n",\
				logfile=logfile,\
				test=test)
		os.system("touch "+outdir+"3_clusters_ok")

	# round 2
	if not os.path.exists(outdir+"4_refine_fasta_ok"):
		write_fasta_files_from_trees.main(fasta=outdir2+"all.fa.cut",\
			treDIR=outdir+"4_refine/",\
			tree_file_ending=".subtree",\
			outDIR=outdir+"4_refine/")
		os.system("touch "+outdir+"4_refine_fasta_ok")
	homodir = outdir+"5_homolog/"
	if not os.path.exists(outdir+"4_refine_ok"):
		refine(curdir=outdir+"4_refine/",nextdir=homodir,\
				prefix="5_homolog",\
				num_cores=num_cores,\
				relative_cut=relative_cut,\
				absolute_cut=absolute_cut,\
				branch_len_cutoff=branch_len_cutoff,\
				min_taxa=min_taxa,\
				mask_para="y",\
				logfile=logfile)
		os.system("touch "+outdir+"4_refine_ok")
	
	# bootstrap the homologs
	if not os.path.exists(outdir+"5_homolog_fasta_ok"):
		# tree files that ends with ".subtree"
		write_fasta_files_from_trees.main(fasta=outdir2+"all.fa.cut",\
			treDIR=homodir,\
			tree_file_ending=".subtree",\
			outDIR=homodir)
		"""
		# tree files that ends with ".mm"
		write_fasta_files_from_trees.main(fasta=outdir2+"all.fa.cut",\
			treDIR=homodir,\
			tree_file_ending=".mm",\
			outDIR=homodir)"""
		os.system("touch "+outdir+"5_homolog_fasta_ok")
	if not os.path.exists(outdir+"5_homolog_bootstrap_ok"):
		for i in os.listdir(homodir):
			if i.endswith(".fa"):
				fasta_to_tree.fasta_to_bs_tree(DIR=homodir,\
						fasta=i,\
						num_cores=num_cores,\
						seqtype="dna")
		os.system("touch "+outdir+"5_homolog_bootstrap_ok")
	print("homologs with bootstrap values (.tre) written to "+homodir)
	
	# get orthologs
	orthodir = outdir+"6_ortho/"
	try: os.stat(orthodir)
	except: os.mkdir(orthodir)
	if not os.path.exists(outdir+"ortho_tre_ok"):
		if args.ortho == "121":
			get_121(indir=homodir,\
				tree_file_ending=".tre",\
				min_taxa=min_taxa,\
				outdir=orthodir,\
				min_bootstrap=80.0)
			os.system("touch "+outdir+"ortho_tre_ok")
			with open(logfile,"a") as f:
				f.write("Filter one-to-one orthologs by minimal bootstrap of 80\n")
		else:
			ingroups = 0
			with open(args.inout) as infile:
				for line in infile:
					if line.startswith("IN\t"): ingroups += 1
			RT(homoDIR=homodir,\
				tree_file_eneding =".tre",\
				outDIR=orthodir,\
				MIN_INGROUP_TAXA=get_min_taxa(ingroups,args.max_mis_taxa),\
				taxon_code_file=args.inout)
			os.system("touch "+outdir+"ortho_tre_ok")
			with open(logfile,"a") as f:
				f.write("RT orthologs\n")
	if not os.path.exists(outdir+"ortho_aln_ok"):
		ortho_to_aln(alndir=homodir,\
			tredir=orthodir,\
			outdir=orthodir,\
			ortho_tree_file_ending=".tre")
		os.system("touch "+outdir+"ortho_aln_ok")
		with open(logfile,"a") as f:
			f.write("Write ortholog alignments from homologs\n")
	min_col_occup = 0.5 if min_taxa < 10 else 0.3 
	if not os.path.exists(outdir+"ortho_cln_ok"):
		phyutility_wrapper.main(DIR=orthodir,\
			min_col_occup=min_col_occup,\
			seqtype="dna")
		os.system("touch "+outdir+"ortho_cln_ok")
		with open(logfile,"a") as f:
			f.write("Trim ortholog alignments by "+str(min_col_occup)+"\n")
	
	# tree inference and jackknife
	# concatenate
	if args.ortho == "121":
		outname = "filter300-"+str(min_taxa)
		if not os.path.exists(outdir+outname+".phy"):
			concatenate(clnDIR=orthodir,\
				numofsitesFilter=300,\
				numoftaxaFilter=min_taxa,\
				seqtype="dna",\
				outname=outname)
	else:
		min_ingroup_taxa = get_min_taxa(ingroups,args.max_mis_taxa)
		outname = "filter300-"+str(min_ingroup_taxa)
		if not os.path.exists(outdir+outname+".phy"):
			concatenate(clnDIR=orthodir,\
				numofsitesFilter=300,\
				numoftaxaFilter=min_ingroup_taxa,\
				seqtype="dna",\
				outname=outname)
	# tree
	if not os.path.exists("RAxML_bestTree."+outname):
		cmd = ["raxml",\
			   "-T",str(num_cores),\
			   "-p","1234",\
			   "-m","GTRCAT",\
			   "-q",outname+".model",\
			   "-s",outname+".phy",\
			   "-n",outname]
		print(" ".join(cmd))
		p = subprocess.Popen(cmd,stdout=subprocess.PIPE)
		out = p.communicate()
		assert p.returncode == 0,"Error raxml"+out[0]
		try:
			os.remove("RAxML_info."+outname)
			os.remove("RAxML_log."+outname)
			os.remove("RAxML_parsimonyTree."+outname)
			os.remove("RAxML_result."+outname)
		except: pass # no need to worry about extra intermediate files
	
	# jackknife
	if not os.path.exists("JK5_trees"):
		jk(indir="./",\
		   num_core=num_cores,\
		   resample_num=5,\
		   seqtype="dna",\
		   replicates=200)
	
	# map jackknife support to species tree
	cmd = ["raxml", "-f","b",\
		   "-t", "RAxML_bestTree."+outname,\
		   "-z", "JK5_trees",
		   "-T", str(num_cores),\
		   "-m", "GTRCAT",\
		   "-n", outname+"_JK5"]
	print(" ".join(cmd))
	p = subprocess.Popen(cmd,stdout=subprocess.PIPE)
	out = p.communicate()
	assert p.returncode == 0,"Error raxml"+out[0]
	try:
		os.remove("RAxML_info."+outname+"_JK5")
		os.remove("RAxML_log."+outname+"_JK5")
		os.remove("RAxML_parsimonyTree."+outname+"_JK5")
		os.remove("RAxML_result."+outname+"_JK5")
		os.remove("RAxML_bipartitionsBranchLabels."+outname+"_JK5")
	except: pass # no need to worry about extra intermediate files
	
	# remove intermediate files
	os.system("rm "+outdir+"*_ok")
	os.remove(outdir+"phyutility.log")

if __name__ == "__main__":
	if len(sys.argv) == 1:
		print("Specify input dir, output dir, and number of threads")
		print("-h for help")
		sys.exit()
		
	main()


