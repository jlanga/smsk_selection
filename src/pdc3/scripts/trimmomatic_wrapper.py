
"""
Wrapper for Trimmomatic
Change paths as needed
"""


import sys,os

APPS_HOME = "/usr/local/bin/" # where Trimmomatic is located
TRIMMOMATIC_CMD = 'java -jar '+APPS_HOME+"/Trimmomatic-0.36/trimmomatic-0.36.jar" # Basic command to call Trimmomatic
TruSeq_ADAPTER = os.path.expanduser("~/Desktop/botany_2018/databases/TruSeq_adapters.fa") # User provided adapters

def trimmomatic_pe(pe_fq1,pe_fq2,num_cores,DIR):
	
	assert os.path.exists(TruSeq_ADAPTER),"Cannot fine the adapter file "+TruSeq_ADAPTER
	
	trim_setting = "ILLUMINACLIP:"+TruSeq_ADAPTER+":2:30:10 SLIDINGWINDOW:4:5 LEADING:5 TRAILING:5 MINLEN:25"
		
	if DIR == ".": DIR = os.getcwd()	
	if os.path.isabs(DIR) == False: DIR = os.path.abspath(DIR)
	if DIR[-1] != "/": DIR += "/"

	
	#path_pe_1, file_pe_1 = (os.path.split(pe_fq1)) #splits the path from the file name
	#path_pe_2, file_pe_2 = (os.path.split(pe_fq2))
	
	path_pe_1, file_pe_1 = os.path.split(pe_fq1)
	pe_fq1_name = str(file_pe_1)
	base_name_pe_1 = pe_fq1_name.split( "." )
	path_pe_2, file_pe_2 = os.path.split(pe_fq2)
	pe_fq2_name = str(file_pe_2)
	base_name_pe_2 = pe_fq2_name.split( "." )

	
	#trimmed_pe1 = (os.path.splitext(file_pe_1)[0])+".paired.trim.fq"
	#trimmed_up1 = (os.path.splitext(file_pe_1)[0])+".unpaired.trim.fq"
	#trimmed_pe2 = (os.path.splitext(file_pe_2)[0])+".paired.trim.fq"
	#trimmed_up2 = (os.path.splitext(file_pe_2)[0])+".unpaired.trim.fq"

	if (os.path.splitext(file_pe_1)[-1]) and (os.path.splitext(file_pe_2)[-1]) == ".gz" :	
		trimmed_pe1 = (base_name_pe_1[0])+".paired.trim.fq.gz"
		trimmed_up1 = (base_name_pe_1[0])+".unpaired.trim.fq.gz"
		trimmed_pe2 = (base_name_pe_2[0])+".paired.trim.fq.gz"
		trimmed_up2 = (base_name_pe_2[0])+".unpaired.trim.fq.gz"
		
	else:
		trimmed_pe1 = (base_name_pe_1[0])+".paired.trim.fq"
		trimmed_up1 = (base_name_pe_1[0])+".unpaired.trim.fq"
		trimmed_pe2 = (base_name_pe_2[0])+".paired.trim.fq"
		trimmed_up2 = (base_name_pe_2[0])+".unpaired.trim.fq"

	
	if os.path.exists(DIR+trimmed_pe1) \
		and os.path.exists(DIR+trimmed_pe2) \
		and os.path.exists(DIR+trimmed_up1) \
		and os.path.exists(DIR+trimmed_up2):
		print(("Found", trimmed_pe1, trimmed_up1, trimmed_pe2, trimmed_up2))
	else:
		cmd = [TRIMMOMATIC_CMD,"PE","-threads",str(num_cores), \
				pe_fq1,pe_fq2, \
				DIR+trimmed_pe1, \
				DIR+trimmed_up1, \
				DIR+trimmed_pe2, \
				DIR+trimmed_up2, \
				trim_setting]
		print((" ".join(cmd)))
		os.system(" ".join(cmd))
	
	assert os.path.exists(DIR+trimmed_pe1) \
			and os.path.exists(DIR+trimmed_pe2) \
			and os.path.exists(DIR+trimmed_up1) \
			and os.path.exists(DIR+trimmed_up2), "Trimmomatic not completed"
			
			
def trimmomatic_se(se_fq,num_cores,DIR):
	
	assert os.path.exists(TruSeq_ADAPTER),"Cannot fine the adapter file "+TruSeq_ADAPTER
	
	trim_setting = "ILLUMINACLIP:"+TruSeq_ADAPTER+":2:30:10 SLIDINGWINDOW:4:5 LEADING:5 TRAILING:5 MINLEN:25"
	
	if DIR == ".": DIR = os.getcwd()	
	if os.path.isabs(DIR) == False: DIR = os.path.abspath(DIR)
	if DIR[-1] != "/": DIR += "/"

	
	#path_se, file_se = (os.path.split(se_fq)) #splits the path from the file name
	#trimmed_se = (os.path.splitext(file_se)[0])+".trim.fq"

	path_se, file_se = (os.path.split(se_fq)) #splits the path from the file name
	se_fq_name = str(file_se)
	base_name_se = se_fq_name.split( "." )

	if (os.path.splitext(file_se)[-1]) == ".gz" :	
		trimmed_se = (base_name_se[0])+".trim.fq.gz"		
	else:
		trimmed_se = (base_name_se[0])+".trim.fq"
		
	if os.path.exists(DIR+trimmed_se):
		print(("Found", trimmed_se))
	else:
		cmd = [TRIMMOMATIC_CMD,"SE","-threads",str(num_cores), \
				se_fq, \
				DIR+trimmed_se, \
				trim_setting]
		print((" ".join(cmd)))
		os.system(" ".join(cmd))
	
	assert os.path.exists(DIR+trimmed_se), "Trimmomatic not completed"
	
if __name__ == "__main__":
	if len(sys.argv) == 4:
		trimmomatic_se(se_fq=sys.argv[1],num_cores=int(sys.argv[2]),DIR=sys.argv[3])
	elif len(sys.argv) == 5:
		trimmomatic_pe(pe_fq1=sys.argv[1],pe_fq2=sys.argv[2],num_cores=int(sys.argv[3]),DIR=sys.argv[4])
	else:
		print ("Usage:")
		print ("For single end reads: python trimmomatic_wrapper.py fastq_se_reads num_cores output_dir ")
		print ("For paired end reads: python trimmomatic_wrapper.py fastq_pe_reads1 fastq_pe_reads2 num_cores output_dir")
		sys.exit(0)

