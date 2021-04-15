import newick3,phylo3,os,sys
from tree_utils import *

def get_121(indir,tree_file_ending,min_taxa,outdir,min_bootstrap=0.0):
	if indir[-1] != "/": indir += "/"
	if outdir[-1] != "/": outdir += "/"
	min_taxa = int(min_taxa)
	min_bootstrap = float(min_bootstrap)
	infile_count, outfile_count = 0,0
	print "Filter one-to-one homologs with average bootstrap of at least",\
		min_bootstrap
	for i in os.listdir(indir):
		if not i.endswith(tree_file_ending): continue
		infile_count += 1
		with open(indir+i,"r") as infile: #only 1 tree in each file
			intree = newick3.parse(infile.readline())
		names = get_front_names(intree)
		num_tips,num_taxa = len(names),len(set(names))
		print "number of tips:",num_tips,"number of taxa:",num_taxa
		if num_tips == num_taxa and num_taxa >= min_taxa:
			if min_bootstrap > 0.0 and not pass_boot_filter(intree,min_bootstrap):
				continue
			print i,"written to out dir"
			outname = i.split(".")[0]+".1to1ortho.tre"
			os.system("cp "+indir+i+" "+outdir+outname)
			outfile_count += 1
	assert infile_count > 0,\
		"No file ends with "+tree_file_ending+" was found in "+indir
	print infile_count,"files read,",outfile_count,"written to",outdir

if __name__ == "__main__":
	if len(sys.argv) != 5:
		print "python filter_1to1_orthologs.py homoTreeDIR tree_file_ending minimal_taxa outDIR"
		sys.exit(0)
	
	indir,tree_file_ending,min_taxa,outdir = sys.argv[1:]
	get_121(indir,tree_file_ending,min_taxa,outdir)

