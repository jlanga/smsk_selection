rule homologs_join_cds:
    input: expand(CDHIT + "{species}.cds", species=SPECIES)
    output: HOMOLOGS + "all.cds"
    conda: "homologs.yml"
    shell:
        """
        cat /dev/null > {output}

        for file in {input}; do
            species=$(basename $file | sed "s/.cds//")
            seqtk seq $file \
            | cut -f 1 -d " " \
            | paste - - \
            | sed "s/^>/>$species@/" \
            | tr "\t" "\n" \
            >> {output}
        done
        """

rule homologs_join_pep:
    input: expand(CDHIT + "{species}.pep", species=SPECIES)
    output: HOMOLOGS + "all.pep"
    conda: "homologs.yml"
    shell:
        """
        cat /dev/null > {output}

        for file in {input}; do
            species=$(basename $file | sed "s/.pep//")
            seqtk seq $file \
            | cut -f 1 -d " " \
            | paste - - \
            | sed "s/^>/>$species@/" \
            | tr "\t" "\n" \
            >> {output}
        done
        """

rule homologs_round1_prepare_msa:
    """
    Transform _ into @, and store it in .aln-cln
    """
    input: ORTHOFINDER + "msa"
    output: touch(HOMOLOGS + "round1_prepare_msa.ok")
    shell: 
        # """
        # mkdir -p {HOMOLOGS_R1}
        # find {input}/ -type f -name "OG*.fa" -exec \
        # bash -c 'sed "s/_/@/" $1 > {HOMOLOGS_R1}/${{1##*/.fa}}.aln-cln' _ {{}} \;
        # """
        """
        mkdir -p {HOMOLOGS_R1}

        find {input}/ -type f -name "OG*.fa" -exec \
        bash -c 'sed "s/_/@/" $1 > {HOMOLOGS_R1}/$(basename $1 .fa).aln-cln' _ {{}} \;
        """

# def aggregate_homologs_round1_prepare(wildcards):
#     checkpoint_trees = checkpoints.orthofinder.get(**wildcards).input[0]
#     trees = expand(
#         ORTHOFINDER + "resolved_gene_trees/{i}.txt",
#         i=glob_wildcards(os.path.join(checkpoint_trees, "{i}.txt")).i
#     )
#     return trees



rule homologs_round1_prepare_trees:
    input: ORTHOFINDER + "gene_trees"
    output: touch(HOMOLOGS + "round1_prepare.ok")
    conda: "homologs.yml"
    shell:
        """
        python src/correct_tree_leaf_names.py {input} .txt {HOMOLOGS_R1} .nwk

        find {HOMOLOGS_R1} -type f -name "OG*_tree.nwk" -exec \
            bash -c 'mv $1 ${{1%_*}}.nwk' _ {{}} \;
        """




rule homologs_round1_treeshrink:
    input:
        HOMOLOGS + "round1_prepare.ok",
        HOMOLOGS + "round1_prepare_msa.ok"
    output: touch(HOMOLOGS + "round1_treeshrink.ok")
    log: HOMOLOGS + "round1_treeshrink.log"
    benchmark: HOMOLOGS + "round1_treeshrink.bmk"
    threads: 32
    params:
        quantile = 0.05
    conda: "homologs.yml"
    shell:
        """
        PATH="bin:$PATH"
        python src/pdc3/scripts/tree_shrink_wrapper.py \
            {HOMOLOGS_R1} \
            .nwk \
            {params.quantile} \
            {HOMOLOGS_R1} \
        2> {log} 1>&2

        find {HOMOLOGS_R1} -name "OG*.ts.tt" -type f -exec \
            bash -c 'mv $1 ${{1%_*}}.ts.tt' _ {{}} \;
        """


rule homologs_round1_mask_tips_by_taxon_id:
    input:
        HOMOLOGS + "round1_treeshrink.ok",
    output: touch(HOMOLOGS + "round1_mask_tips_by_taxon_id.ok")
    params:
        in_dir = HOMOLOGS_R1,
        mask_tips = params["homologs"]["mask_tips"]["mask_paraphyletic"]
    threads: 1
    log: HOMOLOGS + "round1_mask_tips_by_taxon_id.log",
    benchmark: HOMOLOGS + "round1_mask_tips_by_taxon_id.bmk"
    conda: "homologs.yml"
    shell:
        "python src/pdc3/scripts/mask_tips_by_taxonID_transcripts.py "
            "{params.in_dir} "
            "{params.in_dir} "
            "{params.mask_tips} "
        "2> {log} 1>&2"


rule homologs_round1_cut_internal_long_branches:
    input: HOMOLOGS + "round1_mask_tips_by_taxon_id.ok"
    output: touch(HOMOLOGS + "round1_cut_internal_long_branches.ok")
    params:
        in_dir = HOMOLOGS_R1,
        internal_branch_cutoff = params["homologs"]["cut_internal_long_branches"]["internal_branch_cutoff"],
        minimum_taxa = params["homologs"]["cut_internal_long_branches"]["minimum_taxa"]
    threads: 1
    log: HOMOLOGS + "round1_cut_internal_long_branches.log"
    benchmark: HOMOLOGS + "homologs_round1_cut_internal_long_branches.bmk"
    conda: "homologs.yml"
    shell:
        "python src/pdc3/scripts/cut_long_internal_branches.py "
            "{params.in_dir} "
            ".mm "
            "{params.internal_branch_cutoff} "
            "{params.minimum_taxa} "
            "{params.in_dir} "
        "2> {log} 1>&2"


rule homologs_round1_write_fasta_files_from_trees:
    input:
        fasta = HOMOLOGS + "all.pep",
        ok = HOMOLOGS + "round1_cut_internal_long_branches.ok"
    output:
        touch(HOMOLOGS + "round1_write_fasta_files_from_trees.ok")
    params:
        indir = HOMOLOGS_R1,
        outdir = HOMOLOGS_R2
    log: HOMOLOGS + "round1_write_fasta_files_from_trees.log"
    benchmark: HOMOLOGS + "round1_write_fasta_files_from_trees.bmk"
    conda: "homologs.yml"
    shell:
        "python src/pdc3/scripts/write_fasta_files_from_trees.py "
            "{input.fasta} "
            "{params.indir} "
            ".subtree "
            "{params.indir} "
        "2> {log} 1>&2"


rule homologs_round1:
    input:
        rules.homologs_round1_write_fasta_files_from_trees.output


################################################################################
# ROUND 2
################################################################################

rule homologs_round2_prepare:
    input: HOMOLOGS + "round1_write_fasta_files_from_trees.ok"
    output: touch(HOMOLOGS + "round2_prepare.ok")
    shell:
        "mkdir -p {HOMOLOGS_R2}; "
        "ln $(readlink -f {HOMOLOGS_R1}/*rr.fa) {HOMOLOGS_R2}"


rule homologs_round2_fasta_to_tree:
    input: HOMOLOGS + "round2_prepare.ok"
    output:
        touch(HOMOLOGS + "round2_fasta_to_tree.ok")
    params:
        in_dir = HOMOLOGS_R2
    log: HOMOLOGS + "round2_fasta_to_tree.log"
    benchmark: HOMOLOGS + "round2_fasta_to_tree.bmk"
    threads: 32
    conda: "homologs.yml"
    shell:
        """
        PATH="bin:$PATH"
        python src/pdc3/scripts/fasta_to_tree_pxclsq.py \
            {params.in_dir} \
            {threads} \
            aa \
            n \
        2> {log} 1>&2
        """

rule homologs_round2_treeshrink:
    input: HOMOLOGS + "round2_fasta_to_tree.ok"
    output: touch(HOMOLOGS + "round2_treeshrink.ok")
    log: HOMOLOGS + "round2_treeshrink.log"
    benchmark: HOMOLOGS + "round2_treeshrink.bmk"
    threads: 32
    params:
        in_dir = HOMOLOGS_R2,
        quantile = 0.05
    conda: "homologs.yml"
    shell:
        """
        PATH="bin:$PATH"
        python src/pdc3/scripts/tree_shrink_wrapper.py \
            {params.in_dir} \
            .tre \
            {params.quantile} \
            {params.in_dir} \
        2> {log} 1>&2
        """
        # find {HOMOLOGS_R1} -name "OG*.ts.tt" -type f -exec \
        #     bash -c 'mv $1 ${{1%_*} }.ts.tt' _ {{}} \;
        # """


rule homologs_round2_mask_tips_by_taxon_id:
    input: HOMOLOGS + "round2_treeshrink.ok"
    output: touch(HOMOLOGS + "round2_mask_tips_by_taxon_id.ok")
    params:
        in_dir = HOMOLOGS_R2,
        mask_tips = params["homologs"]["mask_tips"]["mask_paraphyletic"]
    threads: 1
    log: HOMOLOGS + "round2_mask_tips_by_taxon_id.log",
    benchmark: HOMOLOGS + "round2_mask_tips_by_taxon_id.bmk"
    conda: "homologs.yml"
    shell:
        "python src/pdc3/scripts/mask_tips_by_taxonID_transcripts.py "
            "{params.in_dir} "
            "{params.in_dir} "
            "{params.mask_tips} "
        "2> {log} 1>&2"


rule homologs_round2_cut_internal_long_branches:
    input: HOMOLOGS + "round2_mask_tips_by_taxon_id.ok"
    output: touch(HOMOLOGS + "round2_cut_internal_long_branches.ok")
    params:
        in_dir = HOMOLOGS_R2,
        internal_branch_cutoff = params["homologs"]["cut_internal_long_branches"]["internal_branch_cutoff"],
        minimum_taxa = params["homologs"]["cut_internal_long_branches"]["minimum_taxa"]
    threads: 1
    log: HOMOLOGS + "round2_cut_internal_long_branches.log"
    benchmark: HOMOLOGS + "homologs_round2_cut_internal_long_branches.bmk"
    conda: "homologs.yml"
    shell:
        "python src/pdc3/scripts/cut_long_internal_branches.py "
            "{params.in_dir} "
            ".mm "
            "{params.internal_branch_cutoff} "
            "{params.minimum_taxa} "
            "{params.in_dir} "
        "2> {log} 1>&2"


checkpoint homologs_round2_write_fasta_files_from_trees:
    input:
        fasta = HOMOLOGS + "all.pep",
        ok = HOMOLOGS + "round2_cut_internal_long_branches.ok"
    output:
        touch(HOMOLOGS + "round2_write_fasta_files_from_trees.ok")
    params:
        indir = HOMOLOGS_R2,
        outdir = HOMOLOGS_R2
    log: HOMOLOGS + "round2_write_fasta_files_from_trees.log"
    benchmark: HOMOLOGS + "round2_write_fasta_files_from_trees.bmk"
    conda: "homologs.yml"
    shell:
        "python src/pdc3/scripts/write_fasta_files_from_trees.py "
            "{input.fasta} "
            "{params.indir} "
            ".subtree "
            "{params.indir} "
        "2> {log} 1>&2"


rule homologs_round2:
    input:
        rules.homologs_round2_write_fasta_files_from_trees.output

def write_inout(samples, output):
    samples\
        .reset_index()\
        [["inout", "species"]]\
        .to_csv(
            path_or_buf=output,
            sep="\t",
            header=False,
            index=False
        )


rule homologs_create_taxa_inout:
    output: HOMOLOGS + "in_out.tsv"
    log: HOMOLOGS + "in_out.log"
    benchmark: HOMOLOGS + "in_out.bmk"
    run:
        samples\
            .reset_index()\
            [["inout", "species"]]\
            .to_csv(
                path_or_buf=output[0],
                sep='\t',
                index=False,
                header=False
            )

rule homologs_prepare_rt:
    input: rules.homologs_round2.input
    output: touch(HOMOLOGS + "prepare_rt.ok")
    params:
        in_dir = HOMOLOGS_R2,
        out_dir = HOMOLOGS_RT
    log: HOMOLOGS + "prepare_rt.log"
    benchmark: HOMOLOGS + "prepare_rt.bmk"
    shell:
        """
        mkdir -p {params.out_dir}
        ln {params.in_dir}/*.tre {params.out_dir}
        """


rule homologs_prune_paralogs_rt:
    input:
        tsv = HOMOLOGS + "in_out.tsv",
        ok = HOMOLOGS + "prepare_rt.ok"
    output:
        ok = touch(HOMOLOGS + "prune_paralogs_rt.ok")
    params:
        in_dir = HOMOLOGS_R2,
        tree_ending = ".tre",
        out_dir = HOMOLOGS_RT,
        minimum_ingroup_taxa = params["homologs"]["root_to_tip"]["minimum_ingroup_taxa"]
    log: HOMOLOGS + "prune_paralogs_rt.log"
    benchmark: HOMOLOGS + "prune_paralogs_rt.bmk"
    conda: "homologs.yml"
    shell:
        """
        python src/pdc3/scripts/prune_paralogs_RT.py \
            {params.in_dir} \
            {params.tree_ending} \
            {params.out_dir} \
            {params.minimum_ingroup_taxa} \
            {input.tsv} \
        2> {log} 1>&2
        """

rule homologs_rt:
    input: HOMOLOGS + "prune_paralogs_rt.ok"



