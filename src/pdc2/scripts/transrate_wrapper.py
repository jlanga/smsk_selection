"Wrapper for transrate"

import sys,os


def transrate_ref(assembly,pe_fq1,pe_fq2,num_cores,DIR,reference):
	
	if DIR == ".": DIR = os.getcwd()	
	if DIR != os.getcwd() : os.chdir(DIR)
	if DIR[-1] != "/": DIR += "/"
	if DIR[0] != "/": DIR = "/" + str(DIR)
	

	path_assembly, file_assembly = (os.path.split(assembly)) #splits the path from the file name
	assembly_base_name = (os.path.splitext(file_assembly)[0])
	results_name = assembly_base_name+"_transrate_results"

	if os.path.exists(DIR+results_name):
	
		print "Found transrate results folder"
	
	else:
		
	
		cmd = ["transrate","--assembly",assembly,
				"--left",pe_fq1, \
				"--right",pe_fq2, \
				"--threads",str(num_cores), \
				"--reference",reference,\
				"--output",results_name]
		print " ".join(cmd)
		os.system(" ".join(cmd))

	assert os.path.exists(DIR+results_name), "Transrate not completed"


def transrate_no_ref(assembly,pe_fq1,pe_fq2,num_cores,DIR):
	
	if DIR == ".": DIR = os.getcwd()	
	if DIR != os.getcwd() : os.chdir(DIR)
	if DIR[-1] != "/": DIR += "/"
	if DIR[0] != "/": DIR = "/" + str(DIR)

	path_assembly, file_assembly = (os.path.split(assembly)) #splits the path from the file name
	assembly_base_name = (os.path.splitext(file_assembly)[0])
	results_name = assembly_base_name+"_transrate_results"


	if os.path.exists(DIR+results_name):
	
		print "Found transrate results folder"
	
	else:


		cmd = ["transrate","--assembly",assembly,
				"--left",pe_fq1, \
				"--right",pe_fq2, \
				"--threads",str(num_cores),\
				"--output",results_name]
		print " ".join(cmd)
		os.system(" ".join(cmd))
		
	print (DIR+results_name)

	assert os.path.exists(DIR+results_name), "Transrate not completed"


if __name__ == "__main__":
	if len(sys.argv) == 6:
		transrate_no_ref(assembly=sys.argv[1],pe_fq1=sys.argv[2],pe_fq2=sys.argv[3],num_cores=int(sys.argv[4]),DIR=sys.argv[5])
	elif len(sys.argv) == 7:
		transrate_ref(assembly=sys.argv[1],pe_fq1=sys.argv[2],pe_fq2=sys.argv[3],num_cores=int(sys.argv[4]),DIR=sys.argv[5],reference=sys.argv[6])
	else:
		print "Usage:"
		print "python transrate_wrapper.py assembly_fasta fastq_filtered_pe_reads1 fastq_filtered_pe_reads2 num_cores output_dir reference(optional) "
		sys.exit(0)



