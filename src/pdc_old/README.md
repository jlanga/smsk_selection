# phylogenomic_dataset_construction


#### This repository contains updated scripts for the [Phylogenomic dataset construction respository](https://bitbucket.org/yangya/phylogenomic_dataset_construction/src/master/) 

The updated scripts cover from the read processing to transcript processing steps. The orthology inference scripts remain the same with only slightly updates.

Citation: Yang, Y. and S.A. Smith. 2014. Orthology inference in non-model organisms using transcriptomes and low-coverage genomes: improving accuracy and matrix occupancy for phylogenomics. Molecular Biology and Evolution. doi: 10.1093/molbev/msu245


##### After cloning the repository you need to change several paths within the next scripts in order to make them work on your local computer:

**extract_sequences.py**: Change the path of CP_DATABASE and MT_DATABASE both files will be in the repository folder called **databases**.

**rcorrector_wrapper.py**: Change the path of APPS_HOME which is where Rcorrector is located on your computer.

**trimmomatic_wrapper.py**: Change the path of APPS_HOME which is where Trimmomatic is located on your computer. Also change TruSeq_ADAPTER, this file will be in the repository folder called **databases**.

**run_chimera_detection.py**: Change the path of SCRIPTS_HOME, this will be the path to the folder **scripts** from the cloned repository.

**transdecoder_wrapper.py**: Change the path of BLASTP_DB_PATH, this will be the path to your custom blast database. One with proteomes of Arabidopsis and Beta is provided in the repository folder called **databases** as db.

## Dependencies needed to run the scripts. __Note that for some dependencies specific versions are needed, so if you use a different version it may not work:__



[Rcorrector](https://github.com/mourisl/Rcorrector)

[Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) Version 0.36 (newer versions should work)

[Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) Version 2.3.3 (newer versions should work)

[FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) Version 0.11.6 (newer versions should work)

[Trinity](https://github.com/trinityrnaseq/trinityrnaseq/wiki) Version 2.5.1 (newer versions may have a conflict with version of Salmon)

[Transrate](http://hibberdlab.com/transrate/) Version 1.0.3 (Problems have been reported with some libraries of Salmon 0.6.2, please check that Salmon works properly. If you have problems with Salmon 0.6.0 you can install Transrate from [here](https://github.com/dfmoralesb/transrate), this will work with Salmon 0.8.2, which needs to be in the path as salmon)

[Corset](https://github.com/Oshlack/Corset) 

[Salmon](https://github.com/COMBINE-lab/salmon/releases) Version v.0.9.1 (This is used for Corset and have not been tested with newer versions). You need to name this version salmon-0.9.1.

[TransDecoder](https://github.com/TransDecoder/TransDecoder/wiki) Version 5.3.0 (older or newer versions probably will affect when shortening the names of translated transcripts)

[BLAST](ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)

[cd-hit](https://github.com/weizhongli/cdhit/releases) Version 4.6.8 (newer versions should work)

[MCL](https://micans.org/mcl/) Version 14-137

[TreeShrink](https://github.com/uym2/TreeShrink) It works now with Version 1.3.2 (older versions won't work)

[RAxML](https://github.com/stamatak/standard-RAxML) Version 8.2.11  (newer versions should work)

[Phyx](https://github.com/FePhyFoFum/phyx)

[mafft](https://mafft.cbrc.jp/alignment/software/) Version 7.307 (newer versions should work)

[Gblocks](http://molevol.cmima.csic.es/castresana/Gblocks.html) Version 0.91b

[FastTree](http://www.microbesonline.org/fasttree/) Version 2.1.10 (newer versions should work)

[Pasta](https://github.com/smirarab/pasta) Version v1.8.2 (newer versions should work)

[Prank](http://wasabiapp.org/software/prank/) Version v.170427

#### **We recommend running [Croco](https://gitlab.mbb.univ-montp2.fr/mbb/CroCo) for each batch of sequence data to check for cross-contamination before running the pipeline.**


### Step 1: Read processing 

**Note that the pipeline is expecting that the reads are names as taxonID_1.fq.gz or taxonID_1.fq for single end reads and taxonID_1.fq.gz and taxonID_2.fq.gz or taxonID_1.fq and taxonID_2.fq for paired end reads. If redas are not name like that the pipeline won't work**

The read processing will do:

1. Random sequencing error correction with [Rcorrector](https://github.com/mourisl/Rcorrector)
2. Removes read pairs that cannot be corrected
3. Remove sequencing adapters and low quality sequences with [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)
4. Filter organelle reads (cpDNA, mtDNA or both) with [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml). Files containing only organelle reads will be produced which can be use to assemble for example the plastomes with [Fast-Plast](https://github.com/mrmckain/Fast-Plast)
5. Runs [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) to check read quality and detect over-represented reads
6. Remove Over-represented sequences

The script **filter_fq.py** will run all this steps for a given mRNA library. For paired end reads:

	python filter_fq.py taxonID_1.fq.gz taxonID_2.fq.gz Order_name genome_to_filter[cp, mt or both] num_cores output_dir
	
The first two arguments are the read files. The Order_name is the plant Order (eg. Caryophyllales) will be used for bowtie2 to create a database to filter the organelle reads and can be replaced with any plant Order (or any taxonomic rank following [NCBI taxonomy](https://www.ncbi.nlm.nih.gov/taxonomy)) where you study group belongs. For a list of available genomes with their correspondence taxonomy check for the [cp_lookout](https://bitbucket.org/yanglab/phylogenomic_dataset_construction/src/bf3cd62fdfa139f7c43cf53815f189e9be264306/databases/chloroplast_NCBI_reference_sequences_OCT_17_2018_lookout.txt?at=master&fileviewer=file-view-default) or [mt_lookout](https://bitbucket.org/yanglab/phylogenomic_dataset_construction/src/bf3cd62fdfa139f7c43cf53815f189e9be264306/databases/mitochondrion_NCBI_reference_sequences_OCT_17_2018_lookout.txt?at=master&fileviewer=file-view-default) tables in the databases folder. For the organelle genome you can especify **cpDNA**, **mtDNA** or **both**. **num_core** is the number of cpus or threads to used. **output_dir** is where all the output files will be saved (any existing directory can be used).
	
For single end reads:
	
	python filter_fq.py taxonID_1.fq.gz Order_name organelle_genome num_cores output_dir
	

To see the argument needed to run a scripts call the script without arguments like:

	python filter_fq.py
	
	Usage:
	For single end reads: python filter_fq.py fastq_se_reads Order_name genome_to_filter[cp, mt or both] num_cores output_dir clean(optional)
	For paired end reads: python filter_fq.py fastq_pe_reads1 fastq_pe_reads2 Order_name genome_to_filter[cp, mt or both] num_cores output_dir clean(optional)

This apply for all scripts described here.

This will produced the filtered filles called **taxonID_1.overep_filtered.fq.gz** and **taxonID_2.overep_filtered.fq.gz** and the files containing the organelle reads called **taxonID_1.org_reads.fq.gz** and **taxonID_1.org_reads.fq.gz**

This script also produces several intermediate files that won't be used anymore and that are pretty large and need to be removed. For this run the same command used previously plus **clean** at the end. This option can be used from the beginning, but I like to make sure that the final output is correct before removing intermediate files. 

	python filter_fq.py taxonID_1.fq.gz taxonID_2.fq.gz Order_name genome_to_filter[cp, mt or both] num_cores output_dir clean


### Step 2: _de novo_ assembly with [Trinity](https://github.com/trinityrnaseq/trinityrnaseq/wiki)
	
Choose a taxonID for each data set. This taxonID will be used throughout the analysis. Use short taxonIDs with 4-6 letters and digits with no special characters. 

The assembly can take several hours (even days, depending of the size of the library) and uses a lot of resources, so make sure to have access to preferably a cluster or a computer with enough computation power. Remember that you need to use the filtered reads for the filtering step.

For paired end read:

	python trinity_wrapper.py taxonID_1.overep_filtered.fq.gz taxonID_2.overep_filtered.fq.gz taxonID num_cores max_memory_GB strand(stranded or non-stranded) output_dir

The first two arguments are the filtered reads. **num_cores** is the number of threads. max_memory_GB is maximun amount of memory RAM (in GB) to used. **strand** specifies if the library is stranded or not (stranded or non-stranded). **output_dir** is the output directory.

For single end reads:

	python trinity_wrapper.py taxonID_1.overep_filtered.fq.gz taxonID num_cores max_memory_GB strand(stranded or non-stranded) output_dir

The output contining the assembled transcriptome will be called **taxonID.Trinity.fasta**


### Step 3: Transcript filtering and Translation using [TransDecoder](https://github.com/TransDecoder/TransDecoder/wiki) with blastp for choosing orfs

Note: The assembly quality and filtering steps using Transrate is only available for paired end libraries because Transrate does not work with single end reads.


#### Assembly quality

Now we're going to perform a _de novo_ assembly quality analysis with [Transrate](http://hibberdlab.com/transrate/)

	python transrate_wrapper.py taxonID.Trinity.fasta taxonID_1.overep_filtered.fq.gz taxonID_2.overep_filtered.fq.gz num_cores output_dir

The results will be found in folder **taxonID.Trinity_transrate_results**. Remember to use the filtered reads for for the assembly quality.


#### Assembly filtering 

We're gonna use the Transrate results to filter "bad" or low supported transcripts, based on three out of four transrate score components. sCnuc ≤0.25; sCcov; ≤0.25; and sCord; ≤0.5. 

For this we will use the **contigs.csv** file found in the transrate output folder from the previous step.

	python filter_transcripts_transrate.py taxonID.Trinity.fasta taxonID.Trinity_transrate_results/taxonID.Trinity/contigs.csv .

This will produced a file containing only the good transcripts and short names (so it doesn't produce a error with blastx in later steps) called taxonID.good_transcripts.short_name.fa 


#### Chimera removal

Now we need to remove chimeric transcripts (using the method from Yang, Y. and S.A. Smith Optimizing de novo assembly of short-read RNA-seq data for phylogenomics. BMC Genomics 2013, 14:328 doi:10.1186/1471-2164-14-328). This a BLAST-based method that depending of the size of the trascriptome can take a up to a couple of hours.

	python run_chimera_detection.py taxonID.good_transcripts.short_name.fa reference_proteome_fasta num_core output_dir
	
This will produced a file called **taxonID.filtered_transcripts.fa** and the taxonID.chimera_transcripts.fa


#### Transcript clustering with [Corset](https://github.com/Oshlack/Corset)

Corset clusters transcripts from the same putative gene based in read share and we use it to extract one representative transcripts per gene. In this case the transcript to be keept is the largest within a cluster.

##### Run corset and extract representative transcript

For paired end reads:

	python corset_wrapper.py taxonID.filtered_transcripts.fa taxonID_1.overep_filtered.fq.gz taxonID_2.overep_filtered.fq.gz num_cores ouput_directory aligner[salmon or bowtie]

Remember to used the filtered reads for this step. For the aligner option **we recommend to use Salmon**. Corset seems to have problems parsing the output of version of that we tested bowtie.

For single end reads:

	python corset_wrapper.py taxonID.filtered_transcripts.fa taxonID_1.overep_filtered.fq.gz num_cores ouput_directory aligner[salmon or bowtie]

This will produced the file **taxonID_salmon-clusters.txt** that is gonna be used to select the representative trascript

	python filter_corset_output.py taxonID.filtered_transcripts.fa taxonID_salmon-clusters.txt output_dir

This will produced the file **taxonID.largest_cluster_transcripts.fa**. This is the final filtered transcript file that will be used for translation.

### Translation with Transdecoder

TransDecoder provides two options for choosing the orfs to retain after scoring candidate orfs. Ya benchmarked blastp vs. PfamAB vs. both and found that since we have closely related, high quality proteomes available for the Caryophyllales, blastp using a custom blast database is much faster (hours vs. days) and more sensitive than Pfam. You should try both options and see the difference for your data sets. Here we use the blastp-only option.

First we need to create a custom blast database using closely-related, for example we have used high quality proteomes from Arabidopsis thaliana and Beta vulgaris:

This can be created using the following commands.

	cat Atha.fa Beta.fa > db
	makeblastdb -in ./db -parse_seqids -dbtype prot -out db
	
Then we need to find candidate orfs, BLAST all candidate orfs against reference we created and output the final orfs preferentially retaining the orfs with blast hits. This will be done with the script transdecoder_wrapper.py

	python transdecoder_wrapper.py taxonID.largest_cluster_transcripts.fa num_cores strand(stranded or non-stranded) output_dir

This will produced the CDS and PEP files **taxonID.cds.fa** and **taxonID.pep.fa**. This are the files that have to be used for homology search.

Notice that the format of the sequences if taxonID@seqID. The special character "@" is used to separate taxonID and seqID.


## Step 4: Clustering

The input sequences for homology inference can be cds or peptides depending on how closely-related the focal taxa are. 

Before homology search, it always worth spending some time to make sure that the sequence names are formated correctly and peptides and cds have matching sequence names. Check for duplicated names, special characters other than digits, letters and "_", all names follow the format taxonID@seqID, and file names are the taxonID. It's good to check especially when some of the data sets were obtained from elsewhere. Most peptides and CDS files from genome annotation contain long names, spaces and special characters and should be eliminated before clustering.

	python check_names.py DIR file_ending
	
	**DIR** is the folder name where the output of Transdecoder is located and **file_ending** is the file extension (eg. fa)

Reduce redundancy. For amino acids:

	cd-hit -i taxonID.fa -o taxonID.fa.cdhit -c 0.995 -n 5 -T <num_cores>

Or alternatively, for cds (use -r 0 since these should all be positive strand after translation):

	cd-hit-est -i taxonID.fa.cds -o code.fa.cds.cdhitest -c 0.99 -n 10 -r 0 -T <num_cores>

All-by-all blast. Copy all the taxonID.fa.cdhit files (or .cdhitest files) into a new directory. Note that it is important to set the maximum number of hits very high (e.g. 1000) to allow the inclusion of all closely related ingroup and outgroup sequences. I usually use an evalue cutoff of 10 so that I don't need to re-run the all-by-all blast again.

Since blastp takes much longer to complete than blastn, I prefer using seperate input fasta files for all-by-all blastp to keep track of progress. I also carry out the makeblastdb step locally before doing blast on a cluster. This is because makeblastdb will check formatting of sequences, duplicate sequence names and special characters in seqIDs etc. It is easier to fix these locally before moving to a cluster.
	
	python all-by-all_blastp.py <DIR> <file_ending> <num_cores>
	cat *.rawblastp >all.rawblast

Or alternatively, if CDS works better when the taxa diverged recently:
	
	cat *.cdhitest >all.fa
	makeblastdb -in all.fa -parse_seqids -dbtype nucl -out all.fa
	blastn -db all.fa -query all.fa -evalue 10 -num_threads <num_cores> -max_target_seqs 1000 -out all.rawblast -outfmt '6 qseqid qlen sseqid slen frames pident nident length mismatch gapopen qstart qend sstart send evalue bitscore'

[Optional] Remove ends of sequences that are not covered by any blast hits from other taxa. Skip this if wish not to cut ends that are fast-evolving, or using sequences from genome annotation.
	
	python cut_seq_ends.py all.fa all.rawblast

Filter raw blast output by hit fraction and prepare input file for mcl. The input can be rawblastp or rawblastn results. I usually use 0.3 or 0.4 for hit_fraction_cutoff when using sequences assembled from RNA-seq depending on how divergent the sequences are. A low hit-fraction cutoff will output clusters with more incomplete sequences and much larger and sparser alignments, whereas a high hit-fraction cutoff gives tighter clusters but ignores incomplete or divergent sequences. For genome data I use a hit_fraction cutoff of 0.5. You can also set IGNORE_INTRASPECIFIC_HITS to be True to avoid recent gene duplications or isoforms forming tight clusters and break off. 

	python blast_to_mcl.py all.rawblast <hit_fraction_cutoff>

The output file is used as input for mcl. Try a few different hit fraction cutoffs and inflation values. My experience is that MCL is robust to minusLogEvalue cutoffs and using a cutoff of 0 or 5 works ok in all the cases I tested. You loose entire clusters of short genes by using a high minusLogEvalue cutoff. Use the smallest inflation value and hit-fraction cutoff value combination that gives alignable clusters. "--te" specifies number of threads, "-I" specifies the inflation value, and -tf 'gq()' specifies minimal -log transformed evalue to consider, and "-abc" specifies the input file format. Here are some example mcl command lines:

	mcl all.rawblast.hit-frac0.4.minusLogEvalue --abc -te 5 -tf 'gq(5)' -I 1.4 -o hit-frac0.4_I1.4_e5
	mcl all.rawblast.hit-frac0.4.minusLogEvalue --abc -te 5 -tf 'gq(5)' -I 2 -o hit-frac0.4_I2_e5
	mcl all.rawblast.hit-frac0.3.minusLogEvalue --abc -te 5 -tf 'gq(5)' -I 1.4 -o hit-frac0.3_I1.4_e5
	mcl all.rawblast.hit-frac0.3.minusLogEvalue --abc -te 5 -tf 'gq(5)' -I 2 -o hit-frac0.3_I2_e5

Write fasta files for each cluster from mcl output. Make a new directory to put the thousands of output fasta files.
	
	mkdir <outDIR>
	python write_fasta_files_from_mcl.py <fasta files with or without ends cut> <mcl_outfile> <minimal_taxa> <outDIR>

Now we have a new directory with fasta files that look like cluster1.fa, cluster2.fa and so on.


## Step 5: Build homolog trees

Make sure that raxml, fasttree, phyx(pxclsq,pxcat,pxs2phy,pxs2nex), TreeShrink, mafft and pasta are properly installed and excutables are in the path, and the executables are named exactly as raxml, fasttree, phyx(pxclsq,pxcat,pxs2phy,pxs2nex) mafft and pasta respectively.

Align each cluster, trim alignment, and infer a tree. Make sure that you have raxml 8 or later installed so that it reads fasta file. With an ealier version of raxml you will get an error message "Problem reading number of species and sites".

For clusters that have less than 1000 sequences, it will be aligned with mafft (--genafpair --maxiterate 1000), trimmed by a minimal column occupancy of 0.1 and tree inference using raxml. For larger clusters it will be aligned with pasta, trimmed by a minimal column occupancy of 0.01 and tree inference using fasttree. The ouput tree files look like clusterID.raxml.tre or clusterID.fasttree.tre for clusters with 1000 or more sequences. 

	python fasta_to_tree_pxclsq.py <fasta dir> <number_cores> dna/aa bootstrap(y/n)

You can visualize some of the trees and alignments. You can see that tips that are 0.4 or more are pretty much junk. There are also some tips that are much longer than near-by tips that are probably results of assembly artifacts.

Trim these spurious tips with [TreeShrink](https://github.com/uym2/TreeShrink)

	python tree_shrink_wrapper.py DIR tree_file_ending quantile

It outputs the tips that were trimmed in the file .txt and the trimmed trees in the files .tt. You would need to test different quantiles to see which one fit better you data. The TreeShrink uses 0.05 as default but this might be too high for some dataset. This produces that when the outgroups have long branches these get cut, so make sure that you check for the output trees and txt files to see wich branches are getting cut and choose a quantile value for you data (although a single quantile value would not work for all trees)

An alternative option is to trim the tips using relative and absolute length cutoffs. For examples, trim tips that are longer than a relative length cutoff and more than 10 times longer than its sister. Also trim tips that are longer than an absolute value. The output tree also ends with ".tt". Keep input and output trees in the same directory.

python trim_tips.py input_tree_dir tree_file_ending relative_cutoff absolute_cutoff

Mask both mono- and (optional) paraphyletic tips that belong to the same taxon. Keep the tip that has the most un-ambiguous charactors in the trimmed alignment. Keep input and output trees in the same directory.

	python mask_tips_by_taxonID_transcripts.py <.tt dir> <aln-cln dir> mask_paraphyletic(y/n)

For phylogenomic data sets that are from annotated genomes, I would only mask monophyletic tips, and keep the sequence with the shortest terminal branch length. Keep input and output trees in the same directory.

	python mask_tips_by_taxonID_genomes.py <.tt dir>

Cut deep paralogs. If interested in building phylogeny a lower (more stringent) long_internal_branch_cutoff should be used. Use a higher (more relaxed) cutoff if interested in homologs to avoid splitting homologs. This works very well with CDS and less effective amino acid squences. For CDS the branch lengths are mostly determined by synonymous distance and are more consistant than for amino acids. Make sure that the indir and outdir are different directories.

	python cut_long_internal_branches.py <input tree dir> <input tree file ending> <internal_branch_length_cutoff> <minimal no. taxa> <outDIR>

Write fasta files from trees. The imput tree file ending should be .subtree

	python write_fasta_files_from_trees.py all.fa <cut tree dir> <tree_file_ending> <outDIR>

Repeat the alignment tree estimation, trimming, masking and cutting deep paralogs. Can use a set of more stringent cutoffs in the second round. After the final round, write fasta files from trees using tree files that ends with .subtree, and estimate the final homolog trees.

Alternatively one can calculate the synonymous distance and use that to guide cutting. However, since we are only trying to get well-aligned clusters for tree inference, choice of length cutoffs here can be somewhat arbitary. 

From here a number of further analyses can be done with the homologs, such as gene tree discordance and back translate peptide alignment to codons with pal2nal and investigate signature of natural selection.


## Step 6: Paralogy pruning to infer orthologs. Use one of the following:

1to1: only look at homologs that are strictly one-to-one. No cutting is carried out.
	
	python filter_1to1_orthologs.py <homologDIR> <tree_file_ending> <minimal_taxa> <outDIR>

MI: prune by maximum inclusion. The long_tip_cutoff here is typically the same as the value used when trimming tips. Set OUTPUT_1to1_ORTHOLOGS to False if wish only to ouput orthologs that is not 1-to-1, for example, when 1-to-1 orthologs have already been analyzed in previous steps.

	python prune_paralogs_MI.py <homologDIR> <tree_file_ending> <relative_long_tip_cutoff>  <absolute_long_tip_cutoff> <minimal_taxa> <outDIR>

MO: prune by using homologs with monophyletic, non-repeating outgroups, reroot and cut paralog from root to tip. If no outgroup, only use those that do not have duplicated taxa. Change the list of ingroup and outgroup names first. Set OUTPUT_1to1_ORTHOLOGS to False if wish only to ouput orthologs that is not 1-to-1

	python prune_paralogs_MO.py <homologDIR> <tree_file_ending> <minimal_taxa> <outDIR>

RT: prune by extracting ingroup clades and then cut paralogs from root to tip. If no outgroup, only use those that do not have duplicated taxa. Compile a list of ingroup and outgroup taxonID, with each line begin with either "IN" or "OUT", followed by a tab, and then the taxonID.

	python prune_paralogs_RT.py <homologDIR> <tree_file_ending> <outDIR> <minimal_ingroup_taxa> <ingroup and outgroup taxonIDs>

Or alternatively, if the input homolog tree is already rooted:

	python prune_paralogs_from_rooted_trees.py <homoTreeDIR> <tree_file_ending> <minimal_taxa> <outDIR>


## Step 7: Visualize matrix occupancy stats and constructing the supermatrix

	python ortholog_occupancy_stats.py <ortho_treDIR>

Read in and rank number of taxa per ortholog from highest to lowest. Plot the ranked number of taxa per ortholog

	a <- as.numeric(read.table("ortho_stats")[,1])
	a <- sort(a, decreasing=TRUE)
	pdf(file="taxon_occupancy.pdf")
	plot(a, type="l", lwd=3, ylab="Number of Taxa in Each Ortholog")
	dev.off()

Check taxon_stats to see if any taxa have unusally low number of genes in the orthologs. Open the file taxon_occupancy.pdf and decide the MIN_TAXA filter. Write new fasta files from ortholog trees

	python write_ortholog_fasta_files.py <fasta file with all seqs> <ortholog tree DIR> outDIR MIN_TAXA

Align final orthologs. Play with a few alignment methods here. Try prank or fsa in addition to mafft or sate for more accurate alignments. Prank tend to create lots of gaps when it is not sure about the homology so make sure to check the alignment visually.

	python prank_wrapper.py <inDIR> <outDIR> <file ending> DNA/aa

Trim alignment. I usually use 0.3 for MIN_COLUMN_OCCUPANCY

	python pxclsq_wrapper.py <inDIR> <MIN_COLUMN_OCCUPANCY> DNA/aa

Or use Gblocks for trimming alignments if the sequences are very divergent (change FILE_ENDING first):

	python pep_gblocks_wrapper.py <inDIR> <outDIR>

Choose the minimal cleaned alignment length and minimal number of taxa filters for whether to include an ortholog in the supermatrix. Concatenate selected cleaned matrices:

	python concatenate_matrices_phyx.py <aln-clnDIR> <numofsites> <numoftaxa> <outfile>

This will output a list of cleaned orthology alignments that passed the filter, a summary of taxon matrix occupancies to check whether any taxon is under represented, and a concatenated matrix in phylip and nexus format
