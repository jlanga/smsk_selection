# smsk_selection: A Snakemake pipeline to find orthologs and marks of positive selection

## 1. Description

This is a pipeline to (briefly described):

1. Predict proteins from transcriptomes (transdecoder),
2. Find orhogroups with OrthoFinder, and methods from Yang et al.
3. Find patterns of positive selection with FastCodeML.
4. Annotate transcripts with transdecoder / trinotate
5. Assess transcriptome completeness with Busco


![smsk_selection pipeline](rulegraph.svg)

## 2. First steps

0. Install [conda](https://conda.io/miniconda.html)

1. Install `snakemake`:

```sh
conda install --yes snakemake
```

3. Clone this repo. In case of error with SSL certificates, add `-c http.sslVerify=false`

```sh
git clone --recursive https://github.com/jlanga/smsk_orthofinder.git
```

4. Compile the necessary dependencies: `phyx`, `guidance` and `fastcodeml`:
```sh
bash src/compile_deps.sh
```

5. Introduce the paths to your samples in `samples.tsv`.

6. Run the pipeline as is:

```
snakemake --use-conda --jobs
```

or run it inside a Docker container:

```
bash src/docker_run.sh -j 4 
```



## 3. File organization

The hierarchy of the folder is the one described in [A Quick Guide to Organizing Computational Biology Projects](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1000424):

```
smsk_selection
├── data: raw data, downloaded fastas, databases,....
├── README.md
├── Snakefile: Pipeline runner
├── results: processed data.
|   ├── busco: SCOs identified
|   ├── cdhit: clustered transcriptome
|   ├── homologs: clustered orthogroups as in Yang et al.
|   ├── orthofinder: clustered orthogroups by orthofinder
|   ├── selection: alignments and positive selection results
|   ├── transcriptome: links to input transcriptomes
|   ├── transdecoder: predicted CDS
|   ├── tree: ML and bayesian species tree from 4fold degenerate sites
|   └── trinotate: transcriptome annotation
└── src: additional source code, tarballs, snakefiles, etc.
```


## 4. Requirements

To run this pipeline it should be only necessary to have `snakemake` and `conda` / `mamba`. They together are able to download the required packages to run each step.

In case of doubt, the `Dockerfile` contains the list of the required packages to install.

## Bibliography

- [A Quick Guide to Organizing Computational Biology Projects](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1000424)

- [Snakemake—a scalable bioinformatics workflow engine](http://bioinformatics.oxfordjournals.org/content/28/19/2520)

- [OrhoFinder](https://github.com/davidemms/OrthoFinder)

- [TransDecoder](https://github.com/TransDecoder/TransDecoder)

- [BioPython](https://github.com/biopython)

- [MiniConda](https://conda.io/miniconda.html)

- [Diamond](https://github.com/bbuchfink/diamond)

- [Hmmer](http://hmmer.org)
