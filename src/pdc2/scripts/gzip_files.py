"""
compress all the files end with a user-specified file ending
"""

import os,sys
import gzip

if __name__ == "__main__":
	if len(sys.argv) != 3:
		print "python gzip_fastq_files.py DIR file_ending"
		sys.exit(0)
	
	DIR = sys.argv[1]
	FILE_ENDING = sys.argv[2]
	done = []
	for i in os.listdir(DIR): done.append(i)
	print done
	
	for i in os.listdir(DIR):
		if i.endswith(FILE_ENDING) and i+".gz" not in done:
			com = "gzip "+DIR+"/"+i
			#com ="tar -cjf "+i+".tar.bz2 "+i
			print com
			os.system(com)

