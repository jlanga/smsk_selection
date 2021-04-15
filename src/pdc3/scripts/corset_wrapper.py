"""
Runs salmon or bowtie to remap read to transcritome assembly from Trinity and cluster transcripts using Corset
Salmon, bowtie and Corset must be in the path

"""

import sys, os

def salmon_index(transcript,num_cores,DIR):

	if DIR == ".": DIR = os.getcwd()	
	if os.path.isabs(DIR) == False: DIR = os.path.abspath(DIR)
	if DIR[-1] != "/": DIR += "/"

	path_transcript, file_transcript = os.path.split(transcript) #splits the path from the file name
	transcript_name = str(file_transcript)
	salmon_base_name = transcript_name.split( "." )
	
	index_name = str(salmon_base_name[0])+"_salmon_index"

	if os.path.exists(DIR+index_name): 
		print("Found salmon index folder for",transcript_name)
	else:
		cmd= ["salmon-0.9.1","index","-t",transcript,"-i",DIR+index_name,"--type quasi","-p",str(num_cores)]
		print((" ".join(cmd)))
		os.system(" ".join(cmd))
			
	assert os.path.exists(DIR+index_name), "salmon-index did not finish"


def salmon_quantify_pe(transcript,index,pe_fq1,pe_fq2,num_cores,DIR):

	if DIR == ".": DIR = os.getcwd()	
	if os.path.isabs(DIR) == False: DIR = os.path.abspath(DIR)
	if DIR[-1] != "/": DIR += "/"

	path_transcript, file_transcript = os.path.split(transcript) #splits the path from the file name
	transcript_name = str(file_transcript)
	salmon_base_name = transcript_name.split( "." )

	index_name = str(salmon_base_name[0])+"_salmon_index"
	quant_name = str(salmon_base_name[0])+"_salmon_quant"

	if os.path.exists(DIR+quant_name): 
		print("Found salmon quantification folder for",salmon_base_name[0])
	else:
		cmd= ["salmon-0.9.1","quant","-i",DIR+index_name,"--dumpEq","--libType A","-p",str(num_cores), \
				"-1",pe_fq1, "-2",pe_fq2, "-o", DIR+quant_name]
		print((" ".join(cmd)))
		os.system(" ".join(cmd))
			
	assert os.path.exists(DIR+quant_name), "salmon-quant did not finish"


def salmon_quantify_se(transcript,index,se_fq,num_cores,DIR):

	if DIR == ".": DIR = os.getcwd()	
	if os.path.isabs(DIR) == False: DIR = os.path.abspath(DIR)
	if DIR[-1] != "/": DIR += "/"

	path_transcript, file_transcript = os.path.split(transcript) #splits the path from the file name
	transcript_name = str(file_transcript)
	salmon_base_name = transcript_name.split( "." )

	index_name = str(salmon_base_name[0])+"_salmon_index"
	quant_name = str(salmon_base_name[0])+"_salmon_quant"

	if os.path.exists(DIR+quant_name): 
		print("Found salmon quantification folder for",salmon_base_name[0])
	else:
		cmd= ["salmon-0.9.1","quant","-i",DIR+index_name,"--dumpEq","--libType A","-p",str(num_cores), \
				"-r",se_fq,"-o", DIR+quant_name]
		print((" ".join(cmd)))
		os.system(" ".join(cmd))
			
	assert os.path.exists(DIR+quant_name), "salmon-quant did not finish"
	
def corset_salmon_eq_classes(transcript,eq_classes,DIR):


	if DIR == ".": DIR = os.getcwd()	
	if os.path.isabs(DIR) == False: DIR = os.path.abspath(DIR)
	if DIR[-1] != "/": DIR += "/"

	path_transcript, file_transcript = os.path.split(transcript) #splits the path from the file name
	transcript_name = str(file_transcript)
	salmon_base_name = transcript_name.split( "." )		
	
	clusters_name = str(salmon_base_name[0])+"_salmon"+"-clusters.txt"
	counts_name = str(salmon_base_name[0])+"_salmon"+"-counts.txt"


	if os.path.exists(DIR+clusters_name) and os.path.exists(DIR+counts_name) : 
		print("Found corset-salmon files for",salmon_base_name[0])
	else:
		cmd = ["corset","-i salmon_eq_classes",eq_classes,"-m 5","-p",(DIR+str(salmon_base_name[0])+"_salmon")]
		print((" ".join(cmd)))
		os.system(" ".join(cmd))
			
	assert os.path.exists(DIR+clusters_name) and os.path.exists(DIR+counts_name), "Corsert did not finish"
	

def run_pe_salmon(transcript,pe_fq1,pe_fq2,num_cores,DIR):


	if DIR == ".": DIR = os.getcwd()	
	if os.path.isabs(DIR) == False: DIR = os.path.abspath(DIR)
	if DIR[-1] != "/": DIR += "/"

	path_transcript, file_transcript = os.path.split(transcript) #splits the path from the file name
	transcript_name = str(file_transcript)
	salmon_base_name = transcript_name.split( "." )

	index_name = str(salmon_base_name[0])+"_salmon_index"
	quant_name = str(salmon_base_name[0])+"_salmon_quant"


	salmon_index(transcript,num_cores,DIR)
	salmon_quantify_pe(transcript,DIR+index_name,pe_fq1,pe_fq2,num_cores,DIR)
	corset_salmon_eq_classes(transcript,str(DIR+quant_name)+"/aux_info/eq_classes.txt",DIR)
	
def run_se_salmon(transcript,se_fq,num_cores,DIR):


	if DIR == ".": DIR = os.getcwd()	
	if os.path.isabs(DIR) == False: DIR = os.path.abspath(DIR)
	if DIR[-1] != "/": DIR += "/"

	path_transcript, file_transcript = os.path.split(transcript) #splits the path from the file name
	transcript_name = str(file_transcript)
	salmon_base_name = transcript_name.split( "." )

	index_name = str(salmon_base_name[0])+"_salmon_index"
	quant_name = str(salmon_base_name[0])+"_salmon_quant"


	salmon_index(transcript,num_cores,DIR)
	salmon_quantify_se(transcript,DIR+index_name,se_fq,num_cores,DIR)
	corset_salmon_eq_classes(transcript,str(DIR+quant_name)+"/aux_info/eq_classes.txt",DIR)



def bowtie_index(transcript,num_cores,DIR):

	if DIR == ".": DIR = os.getcwd()	
	if os.path.isabs(DIR) == False: DIR = os.path.abspath(DIR)
	if DIR[-1] != "/": DIR += "/"
	
	if os.path.isabs(transcript) == False: transcript = os.path.abspath(transcript)
	path_transcript, file_transcript = os.path.split(transcript) #splits the path from the file name
	if path_transcript[-1] != "/": path_transcript += "/"

	transcript_name = str(file_transcript)
	bowtie_base_name = transcript_name.split( "." )
	index_name = DIR+str(bowtie_base_name[0])
	
	index_1_bt = (bowtie_base_name[0])+".1.ebwt"
	index_2_bt = (bowtie_base_name[0])+".2.ebwt"
	index_3_bt = (bowtie_base_name[0])+".3.ebwt"
	index_4_bt = (bowtie_base_name[0])+".4.ebwt"
	index_rev_1_bt = (bowtie_base_name[0])+".rev.1.ebwt"
	index_rev_2_bt = (bowtie_base_name[0])+".rev.2.ebwt"
	
	if os.path.exists(DIR+index_1_bt) \
		and os.path.exists(DIR+index_2_bt) \
		and os.path.exists(DIR+index_3_bt) \
		and os.path.exists(DIR+index_4_bt) \
		and os.path.exists(DIR+index_rev_1_bt) \
		and os.path.exists(DIR+index_rev_2_bt):
		print("Found bowtie index for",transcript_name) 

	else:
		cmd= ["bowtie-build","-q",transcript,index_name,"--threads",str(num_cores)]
		print((" ".join(cmd)))
		os.system(" ".join(cmd))
			
		assert os.path.exists(DIR+index_1_bt) \
			and os.path.exists(DIR+index_2_bt) \
			and os.path.exists(DIR+index_3_bt) \
			and os.path.exists(DIR+index_4_bt) \
			and os.path.exists(DIR+index_rev_1_bt) \
			and os.path.exists(DIR+index_rev_2_bt), "bowtie-build not completed"


def bowtie_align_pe(transcript,index,pe_fq1,pe_fq2,num_cores,DIR):

	if DIR == ".": DIR = os.getcwd()	
	if os.path.isabs(DIR) == False: DIR = os.path.abspath(DIR)
	if DIR[-1] != "/": DIR += "/"

	if os.path.isabs(transcript) == False: transcript = os.path.abspath(transcript)
	path_transcript, file_transcript = os.path.split(transcript) #splits the path from the file name
	if path_transcript[-1] != "/": path_transcript += "/"

	transcript_name = str(file_transcript)
	bowtie_base_name = transcript_name.split( "." )

	sam_file = str((bowtie_base_name[0])+"_bowtie.sam")
	bam_file = str((bowtie_base_name[0])+"_bowtie.bam")
	
	
	if os.path.exists(DIR+bam_file):
		print("Found bowtie bam alignment for",bowtie_base_name[0])

	elif os.path.exists(DIR+sam_file):
		cmd = ["samtools view -S","-b",(DIR+sam_file), ">", (DIR+bam_file)]
		print((" ".join(cmd)))
		os.system(" ".join(cmd))
		
		assert os.path.exists(DIR+bam_file), "Samtools not completed"
				
		if os.path.exists(DIR+bam_file):
			os.remove(DIR+sam_file)

	else:
		cmd= ["bowtie","--all","-k 10","--chunkmbs 2580","--rf --nofw","--threads",str(num_cores), (DIR+str(bowtie_base_name[0])),\
				"-1", pe_fq1, "-2", pe_fq2,"-S",(DIR+sam_file)]
		print((" ".join(cmd)))
		os.system(" ".join(cmd))

		assert os.path.exists(DIR+sam_file), "bowtie-align not completed"

		if os.path.exists(DIR+sam_file):
			cmd = ["samtools view -S","-b",(DIR+sam_file), ">", (DIR+bam_file)]
			print((" ".join(cmd)))
			os.system(" ".join(cmd))
		
		assert os.path.exists(DIR+bam_file), "Samtools not completed"
					
		if os.path.exists(DIR+bam_file):
			os.remove(DIR+sam_file)
			
	
def bowtie_align_se(transcript,index,se_fq,num_cores,DIR):

	if DIR == ".": DIR = os.getcwd()	
	if os.path.isabs(DIR) == False: DIR = os.path.abspath(DIR)
	if DIR[-1] != "/": DIR += "/"

	if os.path.isabs(transcript) == False: transcript = os.path.abspath(transcript)
	path_transcript, file_transcript = os.path.split(transcript) #splits the path from the file name
	if path_transcript[-1] != "/": path_transcript += "/"

	transcript_name = str(file_transcript)
	bowtie_base_name = transcript_name.split( "." )

	sam_file = str((bowtie_base_name[0])+"_bowtie.sam")
	bam_file = str((bowtie_base_name[0])+"_bowtie.bam")
	
	
	if os.path.exists(DIR+bam_file):
			print("Found bowtie bam alignment for",bowtie_base_name[0])

	elif os.path.exists(DIR+sam_file):
		cmd = ["samtools view -S","-b",(DIR+sam_file), ">", (DIR+bam_file)]
		print((" ".join(cmd)))
		os.system(" ".join(cmd))		

		if os.path.exists(DIR+bam_file):
			os.remove(DIR+sam_file)
		
		assert os.path.exists(DIR+bam_file), "Samtools not completed"

	else:
		cmd= ["bowtie","--all","-k 10","--threads",str(num_cores), DIR+str(bowtie_base_name[0]),\
				se_fq, "-S", (DIR+sam_file)]
		print((" ".join(cmd)))
		os.system(" ".join(cmd))
		
		assert os.path.exists(DIR+sam_file), "bowtie-align not completed"

		if os.path.exists(DIR+sam_file):
			cmd = ["samtools view -S","-b",(DIR+sam_file), ">", (DIR+bam_file)]
			print((" ".join(cmd)))
			os.system(" ".join(cmd))

		assert os.path.exists(DIR+bam_file), "Samtools not completed"

		if os.path.exists(DIR+bam_file):
			os.remove(DIR+sam_file)
			


def corset_bowtie(transcript,bam,DIR):


	if DIR == ".": DIR = os.getcwd()	
	if os.path.isabs(DIR) == False: DIR = os.path.abspath(DIR)
	if DIR[-1] != "/": DIR += "/"

	path_transcript, file_transcript = os.path.split(transcript) #splits the path from the file name
	transcript_name = str(file_transcript)
	bowtie_base_name = transcript_name.split( "." )		
	
	clusters_name = str(bowtie_base_name[0])+"_bowtie"+"-clusters.txt"
	counts_name = str(bowtie_base_name[0])+"_bowtie"+"-counts.txt"


	if os.path.exists(DIR+clusters_name) and os.path.exists(DIR+counts_name) : 
		print("Corset files found for",bowtie_base_name[0])
	else:
		cmd = ["corset","-i bam",bam,"-m -5","-p",(DIR+str(bowtie_base_name[0])+"_bowtie")]
		print((" ".join(cmd)))
		os.system(" ".join(cmd))
			
	assert os.path.exists(DIR+clusters_name) and os.path.exists(DIR+counts_name), "Corsert did not finish"


def run_pe_bowtie(transcript,pe_fq1,pe_fq2,num_cores,DIR):


	if DIR == ".": DIR = os.getcwd()	
	if os.path.isabs(DIR) == False: DIR = os.path.abspath(DIR)
	if DIR[-1] != "/": DIR += "/"

	path_transcript, file_transcript = os.path.split(transcript) #splits the path from the file name
	transcript_name = str(file_transcript)
	bowtie_base_name = transcript_name.split( "." )

	index_name = DIR+str(bowtie_base_name[0])
	bam_file = str((bowtie_base_name[0])+"_bowtie.bam")

	bowtie_index(transcript,num_cores,DIR)
	bowtie_align_pe(transcript,index_name,pe_fq1,pe_fq2,num_cores,DIR)
	corset_bowtie(transcript,DIR+bam_file,DIR)
	
	
def run_se_bowtie(transcript,se_fq,num_cores,DIR):


	if DIR == ".": DIR = os.getcwd()	
	if os.path.isabs(DIR) == False: DIR = os.path.abspath(DIR)
	if DIR[-1] != "/": DIR += "/"

	path_transcript, file_transcript = os.path.split(transcript) #splits the path from the file name
	transcript_name = str(file_transcript)
	bowtie_base_name = transcript_name.split( "." )

	index_name = DIR+str(bowtie_base_name[0])
	bam_file = str((bowtie_base_name[0])+"_bowtie.bam")


	bowtie_index(transcript,num_cores,DIR)
	bowtie_align_se(transcript,index_name,se_fq,num_cores,DIR)
	corset_bowtie(transcript,DIR+bam_file,DIR)



if __name__ == "__main__":
	if len(sys.argv) == 6 and sys.argv[5]=="salmon":
		run_se_salmon(transcript=sys.argv[1],se_fq=sys.argv[2],num_cores=int(sys.argv[3]),DIR=sys.argv[4])
	elif len(sys.argv) == 7 and sys.argv[6]=="salmon":
		run_pe_salmon(transcript=sys.argv[1],pe_fq1=sys.argv[2],pe_fq2=sys.argv[3],num_cores=int(sys.argv[4]),DIR=sys.argv[5])
	elif len(sys.argv) == 6 and sys.argv[5]=="bowtie": 	
		run_se_bowtie(transcript=sys.argv[1],se_fq=sys.argv[2],num_cores=int(sys.argv[3]),DIR=sys.argv[4])
	elif len(sys.argv) == 7 and sys.argv[6]=="bowtie":
		run_pe_bowtie(transcript=sys.argv[1],pe_fq1=sys.argv[2],pe_fq2=sys.argv[3],num_cores=int(sys.argv[4]),DIR=sys.argv[5])
	else:
		print ("Usage:")
		print ("For single end reads: python corset_wrapper.py transcript fastq_se_reads num_cores output_dir aligner[salmon or bowtie]")
		print ("For paired end reads: python corset_wrapper.py transcript fastq_pe_reads1 fastq_pe_reads2 num_cores output_dir aligner[salmon or bowtie]")
		sys.exit(0)


	
	



