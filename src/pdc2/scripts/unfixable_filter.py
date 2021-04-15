"""
Remove reads that are flaged as unfixable after Rcorrect 
and removes flag 'cor' from sequence header. 
It writes a log with the number of corrected and unfix reads that are keept and removed.
It take single or paired reads.
It can handle .gz files
It has code from  FilterUncorrectabledPEfastq.py from adam h freedman
"""
import sys, os
import gzip       
from itertools import izip, izip_longest


def grouper(iterable, n, fillvalue=None): #Break data into fixed-length chunks or blocks to parse fq in 4 lines at a time
    args = [iter(iterable)] * n
    return izip_longest(fillvalue=fillvalue, *args) 

def filter_unfix_pe(pe_fq1,pe_fq2,DIR):

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

	
	#if (os.path.splitext(file_pe_1)[-1]) == ".gz" :	
	#	filtered_1 = (os.path.splitext(file_pe_1)[0])+".filtered.fq.gz"
	#	filtered_2 = (os.path.splitext(file_pe_2)[0])+".filtered.fq.gz"
	#else:
	#	filtered_1 = (os.path.splitext(file_pe_1)[0])+".filtered.fq"
	#	filtered_2 = (os.path.splitext(file_pe_2)[0])+".filtered.fq"
	
	if (os.path.splitext(file_pe_1)[-1]) and (os.path.splitext(file_pe_2)[-1]) == ".gz" :	
		filtered_1 = (base_name_pe_1[0])+".fix.fq.gz"
		filtered_2 = (base_name_pe_2[0])+".fix.fq.gz"
	else:
		filtered_1 = (base_name_pe_1[0])+".fix.fq"
		filtered_2 = (base_name_pe_2[0])+".fix.fq"			
	
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

		
		r1_cor_count=0
		r2_cor_count=0
		pair_cor_count=0
		unfix_count=0
		with fqin1 as f1, fqin2 as f2:
		    R1=grouper(f1,4)
		    R2=grouper(f2,4)
		    counter=0
		    for entry in R1:
		        counter+=1
		        if counter%100000==0:
		            print ("%s reads processed" % counter)
		    
		        head1,seq1,placeholder1,qual1=[i.strip() for i in entry]
		        head2,seq2,placeholder2,qual2=[j.strip() for j in R2.next()]
		        
		        if 'unfixable' in head1 or 'unfixable' in head2:
		            unfix_count+=1
		        else:
		            if 'cor' in head1:
		                r1_cor_count+=1
		            if 'cor' in head2:
		                r2_cor_count+=1
		            if 'cor' in head1 or 'cor' in head2:
		                pair_cor_count+=1
		            
		            head1=head1.split('l:')[0][:-1] # keeps label information before the Rcorrector flags (low kmer stat, 'cor' and 'unfixable error')
		            head2=head2.split('l:')[0][:-1] 
		            fq1out.write('%s\n' % '\n'.join([head1,seq1,placeholder1,qual1]))
		            fq2out.write('%s\n' % '\n'.join([head2,seq2,placeholder2,qual2]))
		
		
		log_base_name = pe_fq1_name.split( "_" )
		unfix_log=open((DIR+(log_base_name[0])+'_fix_pe.log'),'w')
		unfix_log.write('total PE reads:%s\nremoved PE reads:%s\nretained PE reads:%s\nR1 corrected:%s\nR2 corrected:%s\npairs corrected:%s\n' % (counter,unfix_count,counter-unfix_count,r1_cor_count,r2_cor_count,pair_cor_count))
		            
		fq1out.close()
		fq2out.close() 

		assert os.path.exists(DIR+filtered_1) \
				and os.path.exists(DIR+filtered_2),"Filtering not completed"

def filter_unfix_se(se_fq,DIR):

	if DIR == ".": DIR = os.getcwd()	
	if os.path.isabs(DIR) == False: DIR = os.path.abspath(DIR)
	if DIR[-1] != "/": DIR += "/"

	
	#path_se, file_se = (os.path.split(se_fq)) #splits the path from the file name
	
	path_se, file_se = (os.path.split(se_fq)) #splits the path from the file name
	se_fq_name = str(file_se)
	base_name_se = se_fq_name.split( "." )


	#if (os.path.splitext(file_se)[-1]) == ".gz" :	
	#	filtered = (os.path.splitext(file_se)[0])+".filtered.fq.gz"
	#else:
	#	filtered = (os.path.splitext(file_se)[0])+".filtered.fq"

	if (os.path.splitext(file_se)[-1]) == ".gz" :	
		filtered = (base_name_se[0])+".fix.fq.gz"
	else:
		filtered = (base_name_se[0])+".fix.fq"

	if os.path.exists(DIR+filtered): 
		print ("Found", filtered)
	else:
		if (os.path.splitext(file_se)[-1]) == ".gz" :
			sein=gzip.open(se_fq, 'rb')
			seout=gzip.open(DIR+filtered,'wb')
		else:
			sein=open(se_fq, 'r')
			seout=open(DIR+filtered,'w')

		se_cor_count=0
		unfix_count=0
		with sein as se:
		    SE=grouper(se,4)
		    counter=0
		    for entry in SE:
		        counter+=1
		        if counter%100000==0:
		            print ("%s reads processed" % counter)
		    
		        head,seq,placeholder,qual=[i.strip() for i in entry]
		        
		        if 'unfixable' in head:
		            unfix_count+=1
		        else:
		            if 'cor' in head:
		                se_cor_count+=1
		
		            head=head.split('l:')[0][:-1] # keeps label information before the Rcorrector flags (low kmer stat, 'cor' and 'unfixable error')
		            seout.write('%s\n' % '\n'.join([head,seq,placeholder,qual]))
		
		log_base_name = se_fq_name.split( "_" )
		unfix_log=open((DIR+(log_base_name[0])+'_fix_se.log'),'w')
		unfix_log.write('total SE reads:%s\nremoved SE reads:%s\nretained SE reads:%s\n' % (counter,unfix_count,counter-unfix_count))
		            
		seout.close()
	assert os.path.exists(DIR+filtered),"Filtering not completed"
 

if __name__ == "__main__":
	if len(sys.argv) == 3:
		filter_unfix_se(se_fq=sys.argv[1],DIR=sys.argv[2])
	elif len(sys.argv) == 4:
		filter_unfix_pe(pe_fq1=sys.argv[1],pe_fq2=sys.argv[2],DIR=sys.argv[3])
	else:
		print ("Usage:")
		print ("For single end reads: python unfixable_filter.py fastq_se_reads output_dir")
		print ("For paired end reads: python unfixable_filter.py fastq_pe_reads1 fastq_pe_reads2 output_dir")
		sys.exit(0)


