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

# rule tree_phyx_trim_cols:
#     input: TREE + "supermatrix.fa"
#     output: TREE + "supermatrix_hq.fa"
#     log: TREE + "supermatrix_hq.log"
#     benchmark: TREE + "supermatrix_hq.bmk"
#     params: params["tree"]["min_occupation"]  # proportion of data that is required to be present
#     shell:
#         """
#         ./bin/pxclsq \
#             --seqf {input} \
#             --outf {output} \
#             --prop {params} \
#         2> {log} 1>&2
#         """



# rule tree_modeltest:
#     input: TREE + "supermatrix.fa"
#     output: TREE + "supe"


# rule tree_raxml:
#     input: TREE + "supermatrix.fa"
#     output: TREE + "tree_raxml.nwk"
#     threads: MAX_THREADS
#     log: TREE + "raxml.log"
#     benchmark: TREE + "raxml.bmk"
#     conda: "tree.yml"
#     shell:
#         """
#         raxmlHPC-PTHREADS \

#         """


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
    input: TREE + "supermatrix.fa"