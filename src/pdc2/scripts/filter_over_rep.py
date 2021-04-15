"""
Removes over-represented sequences detected by FastQC
It has code from  RemoveFastqcOverrepSequenceReads.py from adam h freedman

"""

import sys, os
import gzip
import re
from itertools import izip,izip_longest


def grouper(iterable,fillvalue=None): #Break data into fixed-length chunks or blocks to parse fq in 4 lines at a time
    args = [iter(iterable)] * 4
    return izip_longest(fillvalue=fillvalue, *args) 


def ParseFastqcLog(fastqclog):  #Parse the fastqc_data.txt, find and returns the over-represented sequences
    with open(fastqclog) as fp:
        for result in re.findall('Overrepresented sequences(.*?)END_MODULE', fp.read(), re.S):
            seqs=([i.split('\t')[0] for i in result.split('\n')[2:-1]])
    return seqs     

def seqsmatch(overreplist,read): 
    flag=False
    if overreplist!=[]:
        for seq in overreplist:
            if seq in read:
                flag=True
                break
    return flag

def filter_overep_pe(pe_fq1,pe_fq2,pe_fqc1,pe_fqc2,DIR):

	if DIR == ".": DIR = os.getcwd()	
	if os.path.isabs(DIR) == False: DIR = os.path.abspath(DIR)
	if DIR[-1] != "/": DIR += "/"

	path_pe_1, file_pe_1 = os.path.split(pe_fq1)
	pe_fq1_name = str(file_pe_1)
	base_name_pe_1 = pe_fq1_name.split( "." )
	path_pe_2, file_pe_2 = os.path.split(pe_fq2)
	pe_fq2_name = str(file_pe_2)
	base_name_pe_2 = pe_fq2_name.split( "." )

	if (os.path.splitext(file_pe_1)[-1]) and (os.path.splitext(file_pe_2)[-1]) == ".gz" :	
		filtered_1 = (base_name_pe_1[0])+".overep_filtered.fq.gz"
		filtered_2 = (base_name_pe_2[0])+".overep_filtered.fq.gz"
	else:
		filtered_1 = (base_name_pe_1[0])+".overep_filtered.fq"
		filtered_2 = (base_name_pe_2[0])+".overep_filtered.fq"			
	
	if os.path.exists(DIR+filtered_1) and os.path.exists(DIR+filtered_2): 
		print ("Found", filtered_1, filtered_2)

	else:
		
		if (os.path.splitext(file_pe_1)[-1]) and (os.path.splitext(file_pe_2)[-1]) == ".gz" : 
			fqin1=gzip.open(pe_fq1, 'rb') #finq1 = r1handle
			fqin2=gzip.open(pe_fq2, 'rb')
			fq1out=gzip.open(DIR+filtered_1,'wb') #fq1out = r1out
			fq2out=gzip.open(DIR+filtered_2,'wb')
		else:		
			fqin1=open(pe_fq1, 'r') #finq1 = r1handle
			fqin2=open(pe_fq2, 'r')
			fq1out=open(DIR+filtered_1,'w') #fq1out = r1out
			fq2out=open(DIR+filtered_2,'w')

		leftseqs=ParseFastqcLog(pe_fqc1)
		rightseqs=ParseFastqcLog(pe_fqc2)
		       
		counter=0
		failcounter=0
	
	
		with fqin1 as f1, fqin2 as f2:
			R1=grouper(f1)
			R2=grouper(f2)
			for entry in R1:
				counter+=1
				if counter%100000==0:
					print ("%s reads processed" % counter)
    	
		    	    	head1,seq1,placeholder1,qual1=[i.strip() for i in entry]
    			    	head2,seq2,placeholder2,qual2=[j.strip() for j in R2.next()]
    	    	
    	   		 	flagleft,flagright=seqsmatch(leftseqs,seq1),seqsmatch(rightseqs,seq2)
	    	    	
    	    			if True not in (flagleft,flagright):
    	    				fq1out.write('%s\n' % '\n'.join([head1,seq1,'+',qual1]))
	            			fq2out.write('%s\n' % '\n'.join([head2,seq2,'+',qual2]))
	            		else:
	            			failcounter+=1
	
			print ('total PE reads = %s' % counter)
			print ('retained PE reads = %s' % (counter-failcounter))
			print ('removed PE reads = %s' % failcounter)

		
		log_base_name = pe_fq1_name.split( "_" )
		unfix_log=open((DIR+(log_base_name[0])+'_over_pe.log'),'w')
		unfix_log.write('total PE reads:%s\nremoved PE reads:%s\nretained PE reads:%s\n' % (counter,failcounter,(counter-failcounter)))

		fq1out.close()
		fq2out.close()
		
		assert os.path.exists(DIR+filtered_1) \
			and os.path.exists(DIR+filtered_2),"Over-represented filtering not completed"

		
def filter_overep_se(se_fq,se_fqc,DIR):

	if DIR == ".": DIR = os.getcwd()	
	if os.path.isabs(DIR) == False: DIR = os.path.abspath(DIR)
	if DIR[-1] != "/": DIR += "/"

	path_se, file_se = (os.path.split(se_fq)) #splits the path from the file name
	se_fq_name = str(file_se)
	base_name_se = se_fq_name.split( "." )

	if (os.path.splitext(file_se)[-1]) == ".gz" :	
		filtered = (base_name_se[0])+".overep_filtered.fq.gz"
	else:
		filtered = (base_name_se[0])+".overep_filtered.fq"

	if os.path.exists(DIR+filtered): 
		print ("Found", filtered)

	else:

		if (os.path.splitext(file_se)[-1]) == ".gz" :
			sein=gzip.open(se_fq, 'rb')
			seout=gzip.open(DIR+filtered,'wb')
		else:
			sein=open(se_fq, 'r')
			seout=open(DIR+filtered,'w')

		se_seqs=ParseFastqcLog(se_fqc)
	       
		counter=0
		failcounter=0
	
		with sein as f1:
			R1=grouper(f1)
			for entry in R1:	
		    		counter+=1
		        	if counter%100000==0:
		        		print ("%s reads processed" % counter)
		    	
		    		head1,seq1,placeholder1,qual1=[i.strip() for i in entry]
		        
		        	se_flag=seqsmatch(se_seqs,seq1),
		        
		        	if True not in (se_flag):
		            		seout.write('%s\n' % '\n'.join([head1,seq1,'+',qual1]))
		        	else:
		            		failcounter+=1
	
	
			print ('total SE reads = %s' % counter)
			print ('retained SE reads = %s' % (counter-failcounter))
			print ('removed SE reads = %s' % failcounter)


		log_base_name = se_fq_name.split( "_" )
		unfix_log=open((DIR+(log_base_name[0])+'_over_se.log'),'w')
		unfix_log.write('total SE reads:%s\nremoved SE reads:%s\nretained SE reads:%s\n' % (counter,failcounter,(counter-failcounter)))
		            
		seout.close()

	assert os.path.exists(DIR+filtered),"Over-represented filtering not completed"
 

if __name__ == "__main__":
	if len(sys.argv) == 4:
		filter_overep_se(se_fq=sys.argv[1],se_fqc=sys.argv[2],DIR=sys.argv[3])
	elif len(sys.argv) == 6:
		filter_overep_pe(pe_fq1=sys.argv[1],pe_fq2=sys.argv[2],pe_fqc1=sys.argv[3],pe_fqc2=sys.argv[4],DIR=sys.argv[5])
	else:
		print ("Usage:")
		print ("For single end reads: python over_represented_filter.py fastq_se_reads fastQC_data_se output_dir")
		print ("For paired end reads: python over_represented_filter.py fastq_pe_reads1 fastq_pe_reads2 fastQC_data1 fastQC_data2 output_dir")
		sys.exit(0)


