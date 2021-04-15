"""
Detects and trim abnormally long branches using TreeShrink v1.3.2

TreeShrink must be installed and on path
Phyx (pxrr) must be installed and on path

"""


import sys, os, shutil

def trim(inDIR,tree_file_ending,q,outDIR):

	if os.path.isabs(inDIR) == False: inDIR = os.path.abspath(inDIR)
	if inDIR[-1] != "/": inDIR += "/"
	if os.path.isabs(outDIR) == False: outDIR = os.path.abspath(outDIR)
	if outDIR[-1] != "/": outDIR += "/"
	if outDIR == ".": outDIR = os.getcwd()
	
	
	filecount = 0
	
	#runs treeshrink
	for i in os.listdir(inDIR):
		if i.endswith(tree_file_ending):
			print i
			filecount += 1
			cmd= ["run_treeshrink.py","-t", inDIR+i ,"-c","-m per-gene", "-q "+str(q), "-o",outDIR+i+".ts_dir"]
			print (" ".join(cmd))
			os.system(" ".join(cmd))
			
			
	#go into each TS folder and change extension from 'tree_file_ending' to ts
	
	for f in os.listdir(outDIR):
		if f.endswith(".ts_dir"):
			for t in os.listdir(outDIR+f):
				if t.endswith(tree_file_ending):
					if f[-1] != "/": f += "/"			
					basename = os.path.splitext(outDIR+f+t)[0]
					os.rename(outDIR+f+t, basename + ".ts")
					
					
	
	#moves output files to DIR and delete treeshrink individual folders
	for j in os.listdir(outDIR):
		if j.endswith(".ts_dir"):
			source = outDIR+j
			dest = outDIR
			files = os.listdir(source)
			for f in files:
				shutil.move(source+"/"+f, dest)
			shutil.rmtree(outDIR+j)
	
	
	#removes single quotes from tip labels from treeshrink output trees
	for k in os.listdir(outDIR):
		if k.endswith(".ts"):
			with open(outDIR+k, 'r+') as f:
				content = f.read()
				f.seek(0)
				f.truncate()
				f.write(content.replace("'", ""))
			f.close()
	
	#unroot treeshrink ouput trees
	for l in os.listdir(outDIR):
		if l.endswith(".ts"):
			cmd= ["pxrr","-u","-t", outDIR+l,"-o",outDIR+l+".tt"]
			print (" ".join(cmd))
			os.system(" ".join(cmd))			
			
	#delete ts files
	for m in os.listdir(outDIR):
			if m.endswith(".ts"):
    				os.remove(outDIR+m)
        
            
	assert filecount > 0, \
		"No file end with "+tree_file_ending+" found in "+inDIR
			
			
if __name__ == "__main__":
	if len(sys.argv) != 5:
		print "python tree_shrink_wrapper.py inDIR tree_file_ending quantile outDIR"
		sys.exit(0)

	inDIR,tree_file_ending,q,outDIR = sys.argv[1:]
	trim(inDIR,tree_file_ending,q,outDIR)