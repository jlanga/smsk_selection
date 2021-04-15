"""
From raw fastq to over-represented read filtered fastq files.
This script does error correction with Rcorrector,
filters the unfixable reads after Rcorrector,
remove adapter and low quality sequences with Trimmomatic,
filter organelle reads (cp, mt or both) with Bowtie2
runs FastQC to ckeck quality and detect over-represented reads
and filters over-represented reads

It dependends of indiviudal wrappers that can also be run
separately.
"""

import os,sys,shutil
import rcorrector_wrapper
import unfixable_filter
import trimmomatic_wrapper
import extract_sequences
import filter_organelle_reads
import fastqc_wrapper
import filter_over_rep


def filter_fq_se(se_fq,order,genome,num_cores,DIR):
		
	if DIR == ".": DIR = os.getcwd()	
	if os.path.isabs(DIR) == False: DIR = os.path.abspath(DIR)
	if DIR[-1] != "/": DIR += "/"
	
	path_se, file_se = (os.path.split(se_fq)) #splits the path from the file name
	se_fq_name = str(file_se)
	base_name_se = se_fq_name.split( "." )

	
	#error correction		
	rcorrector_wrapper.rcorrector_se(se_fq,num_cores,DIR)
	
	#filter unfixable reads		
		
	if (os.path.splitext(file_se)[-1]) == ".gz" :	
		corrected_se = (DIR+base_name_se[0])+".cor.fq.gz"
	else:
		corrected_se = (DIR+base_name_se[0])+".cor.fq"
	
	unfixable_filter.filter_unfix_se(corrected_se,DIR)
	
	#run trimmomatic

	if (os.path.splitext(file_se)[-1]) == ".gz" :	
		filtered_se = (DIR+base_name_se[0])+".fix.fq.gz"
	else:
		filtered_se = (DIR+base_name_se[0])+".fix.fq"
	
	trimmomatic_wrapper.trimmomatic_se(filtered_se,num_cores,DIR)
	
	
	#extract organelle sequences
	
	if genome == "cp":
		extract_sequences.extract_order_cp(order,DIR)
	elif genome == "mt":
		extract_sequences.extract_order_mt(order,DIR)
	else: extract_sequences.extract_both_cat(order,DIR)
		
	#filter organelle
	

	if genome == "cp":
		reference = (DIR+order+"_cp.fa")
	elif order == "mt":
		genome == (DIR+order+"_mt.fa")
	else: reference = (DIR+order+"_cp_mt.fa")

	path_se, file_se = (os.path.split(se_fq)) #splits the path from the file name
	se_fq_name = str(file_se)
	base_name_se = se_fq_name.split( "." )

	if (os.path.splitext(file_se)[-1]) == ".gz" :	
		trimmed_se = (DIR+base_name_se[0])+".trim.fq.gz"		
	else:
		trimmed_se = (DIR+base_name_se[0])+".trim.fq"
	
	filter_organelle_reads.build_bt2_index(reference,num_cores,DIR) #build the index
	filter_organelle_reads.filter_organelle_se(reference,trimmed_se,num_cores,DIR)
	
	#run fastQC
	
	if (os.path.splitext(file_se)[-1]) == ".gz" :	
		organelle_filtered_se_reads = (DIR+base_name_se[0])+".org_filtered.fq.gz"
	else:
		organelle_filtered_se_reads = (DIR+base_name_se[0])+".org_filtered.fq"


	fastqc_wrapper.fastQC_se(organelle_filtered_se_reads,num_cores,DIR)
	
	#filter over-represented reads
	
	se_fqc = DIR+(base_name_se[0])+".org_filtered_fastqc/fastqc_data.txt"
	
	filter_over_rep.filter_overep_se(organelle_filtered_se_reads,se_fqc,DIR)


def filter_fq_pe(pe_fq1,pe_fq2,order,genome,num_cores,DIR):	
	
	if DIR == ".": DIR = os.getcwd()	
	if os.path.isabs(DIR) == False: DIR = os.path.abspath(DIR)
	if DIR[-1] != "/": DIR += "/"

	path_pe_1, file_pe_1 = os.path.split(pe_fq1)
	pe_fq1_name = str(file_pe_1)
	base_name_pe_1 = pe_fq1_name.split( "." )
	path_pe_2, file_pe_2 = os.path.split(pe_fq2)
	pe_fq2_name = str(file_pe_2)
	base_name_pe_2 = pe_fq2_name.split( "." )


	#error correction
	
	rcorrector_wrapper.rcorrector_pe(pe_fq1,pe_fq2,num_cores,DIR)
	
	#filter unfixable reads

			
	if (os.path.splitext(file_pe_1)[-1]) and (os.path.splitext(file_pe_2)[-1]) == ".gz" :	
		corrected_1 = (DIR+base_name_pe_1[0])+".cor.fq.gz" 
		corrected_2 = (DIR+base_name_pe_2[0])+".cor.fq.gz" 
	else:
		corrected_1 = (DIR+base_name_pe_1[0])+".cor.fq" #input paired fq1 with extension ".cor.fq"
		corrected_2 = (DIR+base_name_pe_2[0])+".cor.fq" #input paired fq2 with extension ".cor.fq"
	

	unfixable_filter.filter_unfix_pe(corrected_1,corrected_2,DIR)
	
			
	#run trimmomatic
	
		
	if (os.path.splitext(file_pe_1)[-1]) and (os.path.splitext(file_pe_2)[-1]) == ".gz" :	
		filtered_1 = (DIR+base_name_pe_1[0])+".fix.fq.gz"
		filtered_2 = (DIR+base_name_pe_2[0])+".fix.fq.gz"
	else:
		filtered_1 = (DIR+base_name_pe_1[0])+".fix.fq"
		filtered_2 = (DIR+base_name_pe_2[0])+".fix.fq"			
	
		
	trimmomatic_wrapper.trimmomatic_pe(filtered_1,filtered_2,num_cores,DIR)
	
	#extract organelle sequences
	
	if genome == "cp":
		extract_sequences.extract_order_cp(order,DIR)
	elif genome == "mt":
		extract_sequences.extract_order_mt(order,DIR)
	else: extract_sequences.extract_both_cat(order,DIR)
		
	#filter organelle
	
	if genome == "cp":
		reference = (DIR+order+"_cp.fa")
	elif genome == "mt":
		reference = (DIR+order+"_mt.fa")
	else: reference = (DIR+order+"_cp_mt.fa")		
			
	if (os.path.splitext(file_pe_1)[-1]) and (os.path.splitext(file_pe_2)[-1]) == ".gz" :	
		trimmed_pe1 = (DIR+base_name_pe_1[0])+".paired.trim.fq.gz"
		trimmed_pe2 = (DIR+base_name_pe_2[0])+".paired.trim.fq.gz"
		
	else:
		trimmed_pe1 = (DIR+base_name_pe_1[0])+".paired.trim.fq"
		trimmed_pe2 = (DIR+base_name_pe_2[0])+".paired.trim.fq"
		
		
	filter_organelle_reads.build_bt2_index(reference,num_cores,DIR) #build the index
	filter_organelle_reads.filter_organelle_pe(reference,trimmed_pe1,trimmed_pe2,num_cores,DIR)
	
	
	#run fastQC
	
	if (os.path.splitext(file_pe_1)[-1]) and (os.path.splitext(file_pe_2)[-1]) == ".gz" :	
		organelle_filtered_reads_1 = (DIR+base_name_pe_1[0])+".org_filtered.fq.gz"
		organelle_filtered_reads_2 = (DIR+base_name_pe_2[0])+".org_filtered.fq.gz"
	
	else: 
		organelle_filtered_reads_1 = (DIR+base_name_pe_1[0])+".org_filtered.fq"
		organelle_filtered_reads_2 = (DIR+base_name_pe_2[0])+".org_filtered.fq"


	fastqc_wrapper.fastQC_pe(organelle_filtered_reads_1,organelle_filtered_reads_2,num_cores,DIR)
	
	#run over-represented reads

	pe_fqc1 = DIR+(base_name_pe_1[0])+".org_filtered_fastqc/fastqc_data.txt"
	pe_fqc2 = DIR+(base_name_pe_2[0])+".org_filtered_fastqc/fastqc_data.txt"
	
	filter_over_rep.filter_overep_pe(organelle_filtered_reads_1,organelle_filtered_reads_2,pe_fqc1,pe_fqc2,DIR)
	
def clean_se(se_fq,order,genome,DIR):

	print("Cleaning intermediate files...")

	if DIR == ".": DIR = os.getcwd()	
	if os.path.isabs(DIR) == False: DIR = os.path.abspath(DIR)
	if DIR[-1] != "/": DIR += "/"
	
	path_se, file_se = (os.path.split(se_fq)) #splits the path from the file name
	se_fq_name = str(file_se)
	base_name_se = se_fq_name.split( "." )



	if (os.path.splitext(file_se)[-1]) == ".gz" :	
		corrected_se = (base_name_se[0])+".cor.fq.gz"
		filtered_se = (base_name_se[0])+".fix.fq.gz"		
		trimmed_se = (base_name_se[0])+".trim.fq.gz"		
		organelle_filtered_se_reads = (base_name_se[0])+".org_filtered.fq.gz"
	else:
		corrected_se = (base_name_se[0])+".cor.fq"
		filtered_se = (base_name_se[0])+".fix.fq"
		trimmed_se = (base_name_se[0])+".trim.fq"
		organelle_filtered_se_reads = (base_name_se[0])+".org_filtered.fq"

	if genome == "cp":
		reference = (order+"_cp.fa")
	elif genome == "mt":
		reference == (order+"_mt.fa")
	else: reference = (order+"_cp_mt.fa")
	
	
	path_ref, file_ref = (os.path.split(reference)) #splits the path from the file name
	ref_name = str(file_ref)
	base_name_ref = ref_name.split( "." )

	fqc_html_se = (base_name_se[0])+".org_filtered_fastqc.html" 
	fqc_folder_se = (base_name_se[0])+".org_filtered_fastqc"
	fqc_zip_se = (base_name_se[0])+".org_filtered_fastqc.zip" 
	sam_file = (base_name_se[0])+".sam"
	index_1_bt2 = (base_name_ref[0])+".1.bt2"
	index_2_bt2 = (base_name_ref[0])+".2.bt2"
	index_3_bt2 = (base_name_ref[0])+".3.bt2"
	index_4_bt2 = (base_name_ref[0])+".4.bt2"
	index_rev_1_bt2 = (base_name_ref[0])+".rev.1.bt2"
	index_rev_2_bt2 = (base_name_ref[0])+".rev.2.bt2"
	
	log_base_name = se_fq_name.split( "_" )
	unfix_log = (log_base_name[0])+'_fix_se.log'
	over_rep_log =(log_base_name[0])+'_over_se.log'


	os.remove(DIR+corrected_se) 
	os.remove(DIR+filtered_se)
	os.remove(DIR+trimmed_se)
	os.remove(DIR+organelle_filtered_se_reads)
	os.remove(DIR+reference)
	os.remove(DIR+fqc_html_se)  
	shutil.rmtree(DIR+fqc_folder_se)
	os.remove(DIR+fqc_zip_se)
	os.remove(DIR+sam_file)
	os.remove(DIR+index_1_bt2)
	os.remove(DIR+index_2_bt2)
	os.remove(DIR+index_3_bt2)
	os.remove(DIR+index_4_bt2)
	os.remove(DIR+index_rev_1_bt2)
	os.remove(DIR+index_rev_2_bt2)
	os.remove(DIR+unfix_log)
	os.remove(DIR+over_rep_log)
	
	
def clean_pe(pe_fq1,pe_fq2,order,genome,DIR):

	print("Cleaning intermediate files...")	

	if DIR == ".": DIR = os.getcwd()	
	if os.path.isabs(DIR) == False: DIR = os.path.abspath(DIR)
	if DIR[-1] != "/": DIR += "/"

	path_pe_1, file_pe_1 = os.path.split(pe_fq1)
	pe_fq1_name = str(file_pe_1)
	base_name_pe_1 = pe_fq1_name.split( "." )
	path_pe_2, file_pe_2 = os.path.split(pe_fq2)
	pe_fq2_name = str(file_pe_2)
	base_name_pe_2 = pe_fq2_name.split( "." )
	base_name_bowtie = pe_fq1_name.split( "_" )



	if (os.path.splitext(file_pe_1)[-1]) and (os.path.splitext(file_pe_2)[-1]) == ".gz" :	
		corrected_1 = (base_name_pe_1[0])+".cor.fq.gz" 
		corrected_2 = (base_name_pe_2[0])+".cor.fq.gz" 
		filtered_1 = (base_name_pe_1[0])+".fix.fq.gz"
		filtered_2 = (base_name_pe_2[0])+".fix.fq.gz"
		trimmed_pe1 = (base_name_pe_1[0])+".paired.trim.fq.gz"
		trimmed_up1 = (base_name_pe_1[0])+".unpaired.trim.fq.gz"
		trimmed_pe2 = (base_name_pe_2[0])+".paired.trim.fq.gz"
		trimmed_up2 = (base_name_pe_2[0])+".unpaired.trim.fq.gz"
		organelle_filtered_reads_1 = (base_name_pe_1[0])+".org_filtered.fq.gz"
		organelle_filtered_reads_2 = (base_name_pe_2[0])+".org_filtered.fq.gz"
	else:
		corrected_1 = (base_name_pe_1[0])+".cor.fq" #input paired fq1 with extension ".cor.fq"
		corrected_2 = (base_name_pe_2[0])+".cor.fq" #input paired fq2 with extension ".cor.fq"
		filtered_1 = (base_name_pe_1[0])+".fix.fq"
		filtered_2 = (base_name_pe_2[0])+".fix.fq"			
		trimmed_pe1 = (base_name_pe_1[0])+".paired.trim.fq"
		trimmed_up1 = (base_name_pe_1[0])+".unpaired.trim.fq"
		trimmed_pe2 = (base_name_pe_2[0])+".paired.trim.fq"
		trimmed_up2 = (base_name_pe_2[0])+".unpaired.trim.fq"
		organelle_filtered_reads_1 = (base_name_pe_1[0])+".org_filtered.fq"
		organelle_filtered_reads_2 = (base_name_pe_2[0])+".org_filtered.fq"

	if genome == "cp":
		reference = (order+"_cp.fa")
	elif order == "mt":
		genome == (order+"_mt.fa")
	else: reference = (order+"_cp_mt.fa")

	path_ref, file_ref = (os.path.split(reference)) #splits the path from the file name
	ref_name = str(file_ref)
	base_name_ref = ref_name.split( "." )


	fqc_html_1 = (base_name_pe_1[0])+".org_filtered_fastqc.html" 
	fqc_html_2 = (base_name_pe_2[0])+".org_filtered_fastqc.html" 


	fqc_folder_pe1 = (base_name_pe_1[0])+".org_filtered_fastqc"
	fqc_folder_pe2 = (base_name_pe_2[0])+".org_filtered_fastqc"

	fqc_zip_pe1 = (base_name_pe_1[0])+".org_filtered_fastqc.zip" 
	fqc_zip_pe2 = (base_name_pe_2[0])+".org_filtered_fastqc.zip" 

	sam_file = (base_name_bowtie[0])+".sam"
	index_1_bt2 = (base_name_ref[0])+".1.bt2"
	index_2_bt2 = (base_name_ref[0])+".2.bt2"
	index_3_bt2 = (base_name_ref[0])+".3.bt2"
	index_4_bt2 = (base_name_ref[0])+".4.bt2"
	index_rev_1_bt2 = (base_name_ref[0])+".rev.1.bt2"
	index_rev_2_bt2 = (base_name_ref[0])+".rev.2.bt2"
	
	log_base_name = pe_fq1_name.split( "_" )
	unfix_log = (log_base_name[0])+'_fix_pe.log'
	over_rep_log = (log_base_name[0])+'_over_pe.log'

	os.remove(DIR+corrected_1) 
	os.remove(DIR+corrected_2) 
	os.remove(DIR+filtered_1)
	os.remove(DIR+filtered_2)
	os.remove(DIR+trimmed_pe1)
	os.remove(DIR+trimmed_up1)
	os.remove(DIR+trimmed_pe2)
	os.remove(DIR+trimmed_up2)
	os.remove(DIR+organelle_filtered_reads_1)
	os.remove(DIR+organelle_filtered_reads_2)
	os.remove(DIR+reference)
	shutil.rmtree(DIR+fqc_folder_pe1)
	shutil.rmtree(DIR+fqc_folder_pe2)
	os.remove(DIR+fqc_zip_pe1)
	os.remove(DIR+fqc_zip_pe2) 
	os.remove(DIR+fqc_html_1)
	os.remove(DIR+fqc_html_2)
	os.remove(DIR+sam_file)
	os.remove(DIR+index_1_bt2)
	os.remove(DIR+index_2_bt2)
	os.remove(DIR+index_3_bt2)
	os.remove(DIR+index_4_bt2)
	os.remove(DIR+index_rev_1_bt2)
	os.remove(DIR+index_rev_2_bt2)
	os.remove(DIR+unfix_log)
	os.remove(DIR+over_rep_log)
		

if __name__ =="__main__":
	if len(sys.argv) == 6:
		filter_fq_se(se_fq=sys.argv[1],order=sys.argv[2],genome=sys.argv[3],num_cores=int(sys.argv[4]),DIR=sys.argv[5])
	elif len(sys.argv) == 7 and sys.argv[6]== "clean":
		filter_fq_se(se_fq=sys.argv[1],order=sys.argv[2],genome=sys.argv[3],num_cores=int(sys.argv[4]),DIR=sys.argv[5])
		clean_se(se_fq=sys.argv[1],order=sys.argv[2],genome=sys.argv[3],DIR=sys.argv[5])		
	elif len(sys.argv) == 7 and sys.argv[6] != "clean":
		filter_fq_pe(pe_fq1=sys.argv[1],pe_fq2=sys.argv[2],order=sys.argv[3],genome=sys.argv[4],num_cores=int(sys.argv[5]),DIR=sys.argv[6])
	elif len(sys.argv) == 8 and sys.argv[7] == "clean":
		filter_fq_pe(pe_fq1=sys.argv[1],pe_fq2=sys.argv[2],order=sys.argv[3],genome=sys.argv[4],num_cores=int(sys.argv[5]),DIR=sys.argv[6])
		clean_pe(pe_fq1=sys.argv[1],pe_fq2=sys.argv[2],order=sys.argv[3],genome=sys.argv[4],DIR=sys.argv[6])
	else:
		print ("Usage:")
		print ("For single end reads: python filter_fq.py fastq_se_reads Order_name genome_to_filter[cp, mt or both] num_cores output_dir clean(optional)")
		print ("For paired end reads: python filter_fq.py fastq_pe_reads1 fastq_pe_reads2 Order_name genome_to_filter[cp, mt or both] num_cores output_dir clean(optional)")
		sys.exit(0)






