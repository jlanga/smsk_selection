rule tree_prepare:
    input:
        ok = HOMOLOGS + "refine2.ok",
        cds = HOMOLOGS + "all.cds"
    output: touch(TREE + "prepare.ok")
    log: TREE + "prepare.log"
    threads: MAX_THREADS
    benchmark: TREE + "prepare.bmk"
    conda: "tree.yml"
    shell:
        """
        mkdir -p {TREE_ALN}

        ln --symbolic --relative {HOMOLOGS_REFINE2}*.maxalign.cds {TREE_ALN}
        #rename.ul .maxalign.cds .aln-cln {TREE_ALN}/*.maxalign.cds
        """


rule tree_extract_4d_sites:
    input: TREE + "prepare.ok"
    output: touch(TREE + "extract_4d_sistes.ok")
    log: TREE + "extract_4d_sistes.log"
    benchmark: TREE + "extract_4d_sistes.bmk"
    conda: "tree.yml"
    shell:
        """
        parallel \
            --jobs {threads} \
            python2.7 src/extract_4d_sites_cds.py \
                {{.}}.cds \
                {{.}}.4d.cds \
        ::: {TREE_ALN}*.maxalign.cds \
        2> {log}
        """



rule tree_supermatrix:
    input: TREE + "extract_4d_sistes.ok"
    output: TREE + "supermatrix.fa"
    params: TREE + "supermatrix"
    log: TREE + "supermatrix.log"
    benchmark: TREE + "supermatrix.bmk"
    shell:
        """
        PATH="bin:$PATH"

        parallel \
            --jobs {threads} \
            cut -f 1 -d @ "<" {{.}}.cds ">" {{.}}.aln-cln \
        ::: {TREE_ALN}*.4d.cds

        python2.7 src/pdc2/scripts/concatenate_matrices_phyx.py \
            {TREE_ALN} \
            0 \
            3 \
            {params} \
        2> {log} 1>&2
        """

rule tree_phyx_trim_cols:
    input: TREE + "supermatrix.fa"
    output: TREE + "supermatrix_hq.fa"
    log: TREE + "supermatrix_hq.log"
    benchmark: TREE + "supermatrix_hq.bmk"
    params: params["tree"]["min_occupation"]  # proportion of data that is required to be present
    shell:
        """
        ./bin/pxclsq \
            --seqf {input} \
            --outf {output} \
            --prop {params} \
        2> {log} 1>&2
        """



rule tree_modeltest:
    input: TREE + "supermatrix_hq.fa"
    output: TREE + "supermatrix_hq.modeltest-ng.out"
    params: TREE + "supermatrix_hq.modeltest-ng"
    benchmark: TREE + "supermatrix_hq.modeltest-ng.bmk"
    threads: MAX_THREADS
    conda: "tree.yml"
    shell:
        """
        modeltest-ng \
            --datatype nt \
            --input {input} \
            --output {params} \
            --processes {threads} \
            --rngseed 1 \
        2> /dev/null 1>&2
        """


rule tree_raxmlng:
    input: 
        fasta = TREE + "supermatrix_hq.fa",
        model = TREE + "supermatrix_hq.modeltest-ng.out"
    output: TREE + "tree_raxml.nwk"
    log: TREE + "supermatrix_hq.raxmlmg.log"
    benchmark: TREE + "supermatrix_hq.raxmlmg.bmk"
    threads: MAX_THREADS
    params:
        bootstraps = params["tree"]["raxml"]["bootstrap_replicates"],
        prefix = TREE + "supermatrix_hq.raxml"
    conda: "tree.yml"
    shell:
        """
        model=$(grep raxml-ng {input.model} \
            | tail -1 \
            | grep -Eo -- '--model [A-Za-z0-9]+' \
            | cut -f 2 -d " " \
        )

        raxml-ng \
            --msa {input.fasta} \
            --model $model \
            --bootstrap \
            --bs-trees {params.bootstraps} \
            --prefix {params.prefix} \
            --threads {threads} \
        2> {log}
        """


    # """
    # exabayes \
    #     -f {input.phy} \
    #     -s 12345 \
    #     -n exabayes \
    #     .m DNA \
    #     -z \
    #     -w {TREE} \
    #     -C 
    # """

rule tree:
    input: rules.tree_raxmlng.output