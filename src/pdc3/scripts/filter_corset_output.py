"""
It takes the corset cluster.txt output and the transcripts fasta file and
return a file a fasta file with the largest transcript per cluster and a
fasta file with the removed(redundant) sequences
"""

import sys, os
import pandas as pd
from Bio import SeqIO


def filter_corset(transcripts,corset_cluster,DIR):

	if DIR == ".": DIR = os.getcwd()
	if os.path.isabs(DIR) == False: DIR = os.path.abspath(DIR)
	if DIR[-1] != "/": DIR += "/"

	path_transcript, files_transcript = os.path.split(transcripts)
	transcripts_name = str(files_transcript)
	base_name_transcripts = transcripts_name.split( "." )

	largest_cluster_transcripts = base_name_transcripts[0]+".largest_cluster_transcripts.fa"
	redundant_transcripts = base_name_transcripts[0]+".redundant_cluster_transcripts.fa"

	if os.path.exists(DIR+largest_cluster_transcripts) and os.path.exists(DIR+redundant_transcripts):
		print("Largest and reduntant transcript files found for", base_name_transcripts[0])

	else:
		#reads corset cluster files as dataframe(df) and add headers
		clusters_df = pd.read_table(corset_cluster, header=None)
		cluster_columns = "seqid cluster".strip().split(' ')
		clusters_df.columns = cluster_columns
		
		#takes the transcripts file and create a df with sequence name and length
		seqid = []
		length = []
		for rec in SeqIO.parse(transcripts, 'fasta'):
		    name = rec.id
		    seq = rec.seq
		    seqLen = len(rec)
		    seqid.append(name)
		    length.append(seqLen)
		d = {"seqid":seqid,"length":length}
		seq_len_df = pd.DataFrame(d)
		
		#filters previous df to contain only info from the surviving clusters of corset
		seq_len_filtered_df= seq_len_df[seq_len_df['seqid'].isin(clusters_df['seqid'])]
		
		#merge previous df to the corset cluster df so it has cluster sequence name, cluster name and sequence length
		clusters_with_len_df=pd.merge(clusters_df, seq_len_filtered_df, on="seqid", how="left")
		
		#sort previous df by length and drop duplicate clusters(redundant) keeping only the one with the largest sequence 
		largest_cluster_df = clusters_with_len_df.sort_values('length', ascending=False).drop_duplicates('cluster').sort_index()
		
		#returns the df with only the redundant sequences info
		removed_cluster_df = clusters_with_len_df[~clusters_with_len_df['seqid'].isin(largest_cluster_df['seqid'])]
		
		#index transcript sequences
		transcripts = SeqIO.index(transcripts, "fasta")
	
		#writes fasta file and cluster info for the largest sequences from corset clusters
		largest_cluster=[]
		largest_cluster_name = largest_cluster_df["seqid"]
		for i in largest_cluster_name: largest_cluster.append(transcripts[i])
		count = SeqIO.write(largest_cluster,(DIR+largest_cluster_transcripts), "fasta")
		print(("Kept %i largest transcripts from corset clusters" % count))
		largest_cluster_df.to_csv((DIR+base_name_transcripts[0]+".largest_cluster.csv"),index=False)
		
		#writes fasta file and cluster info for redundant sequences from corset clusters
		removed_cluster=[]
		removed_cluster_name =  removed_cluster_df["seqid"]
		for i in removed_cluster_name: removed_cluster.append(transcripts[i])
		count = SeqIO.write(removed_cluster,(DIR+redundant_transcripts), "fasta")
		print(("Removed %i redundant transcripts" % count))
		removed_cluster_df.to_csv((DIR+base_name_transcripts[0]+".redundant_cluster.csv"),index=False)

	assert os.path.exists(DIR+largest_cluster_transcripts) and os.path.exists(DIR+redundant_transcripts), "filter_corset_output did not finish"

if __name__ == "__main__":
	if len(sys.argv) == 4:
		filter_corset(transcripts=sys.argv[1],corset_cluster=sys.argv[2],DIR=sys.argv[3])
	else:
		print ("Usage:")
		print ("python filter_corset_output.py transcripts_fasta corset_cluster_file output_dir")
		sys.exit(0)


