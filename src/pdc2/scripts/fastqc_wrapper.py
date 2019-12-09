"Wrapper for fastQC"

import sys,os


def fastQC_pe(pe_fq1,pe_fq2,num_cores,DIR):

	if DIR == ".": DIR = os.getcwd()	
	if os.path.isabs(DIR) == False: DIR = os.path.abspath(DIR)
	if DIR[-1] != "/": DIR += "/"

	#path_pe_1, file_pe_1 = (os.path.split(pe_fq1)) #splits the path from the file name
	#path_pe_2, file_pe_2 = (os.path.split(pe_fq2))
	#fqc_html_1 = (os.path.splitext(file_pe_1)[0])+".fastqc.html" 
	#fqc_html_2 = (os.path.splitext(file_pe_2)[0])+".fastqc.html" 

	path_pe_1, file_pe_1 = os.path.split(pe_fq1)
	pe_fq1_name = str(file_pe_1)
	base_name_pe_1 = pe_fq1_name.split( "." )
	path_pe_2, file_pe_2 = os.path.split(pe_fq1)	
	pe_fq2_name = str(file_pe_2)
	base_name_pe_2 = pe_fq2_name.split( "." )

	fqc_html_1 = (base_name_pe_1[0])+".org_filtered_fastqc.html" 
	fqc_html_2 = (base_name_pe_2[0])+".org_filtered_fastqc.html" 


	if os.path.exists(DIR+fqc_html_1) and os.path.exists(DIR+fqc_html_1): 
		print ("Found", fqc_html_1, fqc_html_1)
	else:
		cmd = ["fastqc",pe_fq1,pe_fq2,"-t",str(num_cores),"-o",DIR,"--extract"]
		print (" ".join(cmd))
		os.system(" ".join(cmd))
	assert os.path.exists(DIR+fqc_html_1) \
			and os.path.exists(DIR+fqc_html_2),"fastQC did not finished"
			
			
def fastQC_se(se_fq,num_cores,DIR):

	if DIR == ".": DIR = os.getcwd()	
	if os.path.isabs(DIR) == False: DIR = os.path.abspath(DIR)
	if DIR[-1] != "/": DIR += "/"
	

	path_se, file_se = (os.path.split(se_fq)) #splits the path from the file name
	pe_fq_name = str(file_se)
	base_name_se = pe_fq_name.split( "." )
	
	fqc_html_se = (base_name_se[0])+".org_filtered_fastqc.html" 


	if os.path.exists(DIR+fqc_html_se): 
		print ("Found", fqc_html_se)
	else:
		cmd = ["fastqc",se_fq,"-t",str(num_cores),"-o",DIR,"--extract"]
		print (" ".join(cmd))
		os.system(" ".join(cmd))
	assert os.path.exists(DIR+fqc_html_se),"fastQC did not finished"
	

if __name__ == "__main__":
	if len(sys.argv) == 4:
		fastQC_se(se_fq=sys.argv[1],num_cores=int(sys.argv[2]),DIR=sys.argv[3])
	elif len(sys.argv) == 5:
		fastQC_pe(pe_fq1=sys.argv[1],pe_fq2=sys.argv[2],num_cores=int(sys.argv[3]),DIR=sys.argv[4])
	else:
		print ("Usage:")
		print ("For single end reads: python fastqc_wrapper.py fastq_se_reads num_cores output_dir")
		print ("For paired end reads: python fastqc_wrapper.py fastq_pe_reads1 fastq_pe_reads2 num_cores output_dir")
		sys.exit(0)




