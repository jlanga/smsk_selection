"""
Filter organelle sequences given a chloroplast or mitochondrial reference genomes
Dependencies:bowtie2
It can handlles *gz files. 
"""

import sys,os

def build_bt2_index(reference,num_cores,DIR):

	if DIR == ".": DIR = os.getcwd()	
	if os.path.isabs(DIR) == False: DIR = os.path.abspath(DIR)
	if DIR[-1] != "/": DIR += "/"
	if os.path.isabs(reference) == False: reference = os.path.abspath(reference)


	path_ref, file_ref = (os.path.split(reference)) #splits the path from the file name
	ref_name = str(file_ref)
	base_name_ref = ref_name.split( "." )

	
	
	if path_ref[-1] != "/": path_ref += "/"
	
	index_1_bt2 = (base_name_ref[0])+".1.bt2"
	index_2_bt2 = (base_name_ref[0])+".2.bt2"
	index_3_bt2 = (base_name_ref[0])+".3.bt2"
	index_4_bt2 = (base_name_ref[0])+".4.bt2"
	index_rev_1_bt2 = (base_name_ref[0])+".rev.1.bt2"
	index_rev_2_bt2 = (base_name_ref[0])+".rev.2.bt2"
	
	if os.path.exists(DIR+index_1_bt2) \
		and os.path.exists(DIR+index_2_bt2) \
		and os.path.exists(DIR+index_3_bt2) \
		and os.path.exists(DIR+index_4_bt2) \
		and os.path.exists(DIR+index_rev_1_bt2) \
		and os.path.exists(DIR+index_rev_2_bt2):
		print "Found bowtie2 index for",file_ref 
	
	else:
		cmd = ["bowtie2-build",reference,(DIR+(base_name_ref[0])),"--threads",str(num_cores),"-q"]
		print (" ".join(cmd))
		os.system(" ".join(cmd))
					
	assert os.path.exists(DIR+index_1_bt2) \
			and os.path.exists(DIR+index_2_bt2) \
			and os.path.exists(DIR+index_3_bt2) \
			and os.path.exists(DIR+index_4_bt2) \
			and os.path.exists(DIR+index_rev_1_bt2) \
			and os.path.exists(DIR+index_rev_2_bt2), "bowtie2-build not completed"
		
	
def filter_organelle_pe(reference,pe_fq1,pe_fq2,num_cores,DIR):

	if DIR == ".": DIR = os.getcwd()	
	if os.path.isabs(DIR) == False: DIR = os.path.abspath(DIR)
	if DIR[-1] != "/": DIR += "/"

	path_ref, file_ref = (os.path.split(reference)) #splits the path from the file name
	ref_name = str(file_ref)
	base_name_ref = ref_name.split( "." )
	
	path_pe_1, file_pe_1 = os.path.split(pe_fq1)
	pe_fq1_name = str(file_pe_1)
	base_name_assert_1 = pe_fq1_name.split( "." )
	
	path_pe_2, file_pe_2 = os.path.split(pe_fq2)
	pe_fq2_name = str(file_pe_2)
	base_name_assert_2 = pe_fq2_name.split( "." )
	
	base_name_bowtie = pe_fq1_name.split( "_" )
	
	if (os.path.splitext(file_pe_1)[-1]) and (os.path.splitext(file_pe_2)[-1]) == ".gz" :	
	
		organelle_filtered_reads_1 = (base_name_assert_1[0])+".org_filtered.fq.gz"
		organelle_filtered_reads_2 = (base_name_assert_2[0])+".org_filtered.fq.gz"
		organelle_reads_1 = (base_name_assert_1[0])+".org_reads.fq.gz"
		organelle_reads_2 = (base_name_assert_2[0])+".org_reads.fq.gz"


		if os.path.exists(DIR+organelle_filtered_reads_1) \
			and os.path.exists(DIR+organelle_filtered_reads_2) \
			and os.path.exists(DIR+organelle_reads_1) \
			and os.path.exists(DIR+organelle_reads_2):
			print "Found bowtie2 organelle filtering files found for",(base_name_bowtie[0])
	
		else:
			cmd = ["bowtie2","-x",(DIR+(base_name_ref[0])),"-1",pe_fq1,"-2",pe_fq2,"--very-sensitive-local","--threads",str(num_cores), \
				"--un-conc-gz",(DIR+(base_name_bowtie[0])+"_%.org_filtered.fq.gz"),"--al-conc-gz",(DIR+(base_name_bowtie[0])+"_%.org_reads.fq.gz"),"-S",(DIR+(base_name_bowtie[0])+".sam")]
			print (" ".join(cmd))
			os.system(" ".join(cmd))
	
		assert os.path.exists(DIR+organelle_filtered_reads_1) \
				and os.path.exists(DIR+organelle_filtered_reads_2) \
				and os.path.exists(DIR+organelle_reads_1) \
				and os.path.exists(DIR+organelle_reads_2),"bowtie2 not completed"


	else: 

		organelle_filtered_reads_1 = (base_name_assert_1[0])+".org_filtered.fq"
		organelle_filtered_reads_2 = (base_name_assert_2[0])+".org_filtered.fq"
		organelle_reads_1 = (base_name_assert_1[0])+".org_reads.fq"
		organelle_reads_2 = (base_name_assert_2[0])+".org_reads.fq"


		if os.path.exists(DIR+organelle_filtered_reads_1) \
			and os.path.exists(DIR+organelle_filtered_reads_2) \
			and os.path.exists(DIR+organelle_reads_1) \
			and os.path.exists(DIR+organelle_reads_2):
			print "Found bowtie2 organelle filtering files found for",(base_name_bowtie[0])
	
		else:
			cmd = ["bowtie2","-x",(DIR+(base_name_ref[0])),"-1",pe_fq1,"-2",pe_fq2,"--very-sensitive-local","--threads",str(num_cores), \
				"--un-conc",(DIR+(base_name_bowtie[0])+"_%.org_filtered.fq"),"--al-conc",(DIR+(base_name_bowtie[0])+"_%.org_reads.fq"),"-S",(DIR+(base_name_bowtie[0])+".sam")]
			print (" ".join(cmd))
			os.system(" ".join(cmd))
	
		assert os.path.exists(DIR+organelle_filtered_reads_1) \
				and os.path.exists(DIR+organelle_filtered_reads_2) \
				and os.path.exists(DIR+organelle_reads_1) \
				and os.path.exists(DIR+organelle_reads_2),"bowtie2 not completed"


def filter_organelle_se(reference,se_fq,num_cores,DIR):

	if DIR == ".": DIR = os.getcwd()	
	if os.path.isabs(DIR) == False: DIR = os.path.abspath(DIR)
	if DIR[-1] != "/": DIR += "/"
	
	path_ref, file_ref = (os.path.split(reference)) #splits the path from the file name
	ref_name = str(file_ref)
	base_name_ref = ref_name.split( "." )

	path_se, file_se = (os.path.split(se_fq)) #splits the path from the file name
	se_name = str(file_se)
	base_name_assert_se = se_name.split( "." )

	if (os.path.splitext(file_se)[-1]) == ".gz" :	

		organelle_filtered_se_reads = (base_name_assert_se[0])+".org_filtered.fq.gz"
		organelle_reads_se = (base_name_assert_se[0])+".org_reads.fq.gz"

	
		if os.path.exists(DIR+organelle_filtered_se_reads) \
			and os.path.exists(DIR+organelle_reads_se):
			print "Found bowtie2 organelle filtering files found for",(base_name_assert_se[0])
	
		else:
			cmd = ["bowtie2","-x",(DIR+(base_name_ref[0])),"-U",se_fq,"--very-sensitive-local","--threads",str(num_cores), \
				"--un-gz",(DIR+(base_name_assert_se[0])+".org_filtered.fq.gz"),"--al-gz", \
				(DIR+(base_name_assert_se[0]))+".org_reads.fq.gz","-S",(DIR+(base_name_assert_se[0]))+".sam"]
			print (" ".join(cmd))
			os.system(" ".join(cmd))
	
		assert os.path.exists(DIR+organelle_filtered_se_reads) \
				and os.path.exists(DIR+organelle_reads_se),"bowtie2 not completed"


	else:
	
		organelle_filtered_se_reads = (base_name_assert_se[0])+".org_filtered.fq"
		organelle_reads_se = (base_name_assert_se[0])+".org_reads.fq"

	
		if os.path.exists(DIR+organelle_filtered_se_reads) \
			and os.path.exists(DIR+organelle_reads_se):
			print "Found bowtie2 organelle filtering files found for",(base_name_assert_se[0])
	
		else:
			cmd = ["bowtie2","-x",(DIR+(base_name_ref[0])),"-U",se_fq,"--very-sensitive-local","--threads",str(num_cores), \
				"--un",(DIR+(base_name_assert_se[0])+".org_filtered.fq"),"--al", \
				(DIR+(base_name_assert_se[0]))+".org_reads.fq","-S",(DIR+(base_name_assert_se[0]))+".sam"]
			print (" ".join(cmd))
			os.system(" ".join(cmd))
	
		assert os.path.exists(DIR+organelle_filtered_se_reads) \
				and os.path.exists(DIR+organelle_reads_se),"bowtie2 not completed"

	

if __name__ == "__main__":
	if len(sys.argv) == 6:
		build_bt2_index(reference=sys.argv[1],num_cores=sys.argv[4],DIR=sys.argv[5])
		filter_organelle_pe(reference=sys.argv[1],pe_fq1=sys.argv[2],pe_fq2=sys.argv[3],num_cores=int(sys.argv[4]),DIR=sys.argv[5])
	elif len(sys.argv) == 5:
		build_bt2_index(reference=sys.argv[1],num_cores=sys.argv[3],DIR=sys.argv[4])
		filter_organelle_se(reference=sys.argv[1],se_fq=sys.argv[2],num_cores=int(sys.argv[3]),DIR=sys.argv[4])
	else:
		print ("Usage:")
		print ("For paired end reads: python filter_organelle_reads.py organelle_reference_file pe_fq1 pe_fq2 num_cores output_dir")
		print ("For single end reads: python filter_organelle_reads.py organelle_reference_file se_fq num_cores output_dir")
		sys.exit(0)
	

	

