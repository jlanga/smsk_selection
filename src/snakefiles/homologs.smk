INTERNAL_BRANCH_CUTOFF = params["homologs"]["cut_internal_long_branches"]["internal_branch_cutoff"],
MINIMUM_TAXA = params["homologs"]["cut_internal_long_branches"]["minimum_taxa"]
MASK_PARAPHYLETIC_TIPS = params["homologs"]["mask_tips"]["mask_paraphyletic"]
MINIMUM_INGROUP_TAXA = params["homologs"]["root_to_tip"]["minimum_ingroup_taxa"]

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
        """
        mkdir -p {HOMOLOGS_R1}

        find {input}/ -type f -name "OG*.fa" -exec \
        bash -c 'sed "s/_/@/" $1 > {HOMOLOGS_R1}/$(basename $1 .fa).aln-cln' _ {{}} \;
        """



rule homologs_round1_prepare_trees:
    """
    Puts the trees in the input folder into the output one, while correcting the
    leaf names for the scripts in pdc2
    """
    input: ORTHOFINDER + "gene_trees"
    output: touch(HOMOLOGS + "round1_prepare.ok")
    conda: "homologs.yml"
    shell:
        """
        python src/correct_tree_leaf_names.py {input} .txt {HOMOLOGS_R1} .nwk

        find {HOMOLOGS_R1} -type f -name "OG*_tree.nwk" -exec \
            bash -c 'mv $1 ${{1%_*}}.nwk' _ {{}} \;
        """


# rule homologs_round1_treeshrink:
#     """
#     Run treeshrink in every tree.
#     Input: .nwk
#     Output: *.ts.tt
#     Rename from *.ts.tt to 
#     """
#     input:
#         HOMOLOGS + "round1_prepare.ok",
#         HOMOLOGS + "round1_prepare_msa.ok"
#     output: touch(HOMOLOGS + "round1_treeshrink.ok")
#     log: HOMOLOGS + "round1_treeshrink.log"
#     benchmark: HOMOLOGS + "round1_treeshrink.bmk"
#     threads: MAX_THREADS
#     params:
#         quantile = 0.05
#     conda: "homologs.yml"
#     shell:
#         """
#         PATH="bin:$PATH"
#         python src/pdc2/scripts/tree_shrink_wrapper.py \
#             {HOMOLOGS_R1} \
#             .nwk \
#             {params.quantile} \
#             {HOMOLOGS_R1} \
#             {threads} \
#         2> {log} 1>&2
#         rm -rf phyx.logfile
#         """

rule homologs_round1_trim_tips:
    """Trim tips via the old method"""
    input:
        HOMOLOGS + "round1_prepare.ok",
        HOMOLOGS + "round1_prepare_msa.ok"
    output: touch(HOMOLOGS + "round1_trim_tips.ok")
    log: HOMOLOGS + "round1_trim_tips.log"
    benchmark: HOMOLOGS + "round1_trim_tips.bmk"
    params:
        relative_cutoff=params["homologs"]["trim_tips"]["relative_cutoff"],
        absolute_cutoff=params["homologs"]["trim_tips"]["absolute_cutoff"]
    conda: "homologs.yml"
    shell:
        """
        python src/pdc2/scripts/trim_tips.py \
            {HOMOLOGS_R1} \
            .nwk \
            {params.relative_cutoff} \
            {params.absolute_cutoff} \
        2> {log} 1>&2
        """


rule homologs_round1_mask_tips_by_taxon_id:
    """
    Result:
    """    
    input: HOMOLOGS + "round1_trim_tips.ok",
    output: touch(HOMOLOGS + "round1_mask_tips_by_taxon_id.ok")
    params:
        in_dir = HOMOLOGS_R1,
        mask_tips = MASK_PARAPHYLETIC_TIPS
    threads: 1
    log: HOMOLOGS + "round1_mask_tips_by_taxon_id.log",
    benchmark: HOMOLOGS + "round1_mask_tips_by_taxon_id.bmk"
    conda: "homologs.yml"
    shell:
        "python src/pdc2/scripts/mask_tips_by_taxonID_transcripts.py "
            "{params.in_dir} "
            "{params.in_dir} "
            "{params.mask_tips} "
        "2> {log} 1>&2"


rule homologs_round1_cut_internal_long_branches:
    """
    Result: OG\d+_\d
    """
    input: HOMOLOGS + "round1_mask_tips_by_taxon_id.ok"
    output: touch(HOMOLOGS + "round1_cut_internal_long_branches.ok")
    params:
        in_dir = HOMOLOGS_R1,
        internal_branch_cutoff = INTERNAL_BRANCH_CUTOFF,
        minimum_taxa = MINIMUM_TAXA
    threads: 1
    log: HOMOLOGS + "round1_cut_internal_long_branches.log"
    benchmark: HOMOLOGS + "round1_cut_internal_long_branches.bmk"
    conda: "homologs.yml"
    shell:
        "python src/pdc2/scripts/cut_long_internal_branches.py "
            "{params.in_dir} "
            ".mm "
            "{params.internal_branch_cutoff} "
            "{params.minimum_taxa} "
            "{params.in_dir} "
        "2> {log} 1>&2"


rule homologs_round1_write_fasta_files_from_trees:
    """
    results: OG\đ+_\d+rr.fa
    """
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
        """
        python src/pdc2/scripts/write_fasta_files_from_trees.py \
            {input.fasta} \
            {params.indir} \
            .subtree \
            {params.indir} \
        2> {log} 1>&2

        find {HOMOLOGS_R1} -name "OG*rr.fa" -type f -exec \
            bash -c 'mv $1 ${{1%rr.fa}}.fa' _ {{}} \;
        """


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
        "ln $(readlink -f {HOMOLOGS_R1}/*.fa) {HOMOLOGS_R2}"


rule homologs_round2_fasta_to_tree:
    input: HOMOLOGS + "round2_prepare.ok"
    output:
        touch(HOMOLOGS + "round2_fasta_to_tree.ok")
    params:
        in_dir = HOMOLOGS_R2
    log: HOMOLOGS + "round2_fasta_to_tree.log"
    benchmark: HOMOLOGS + "round2_fasta_to_tree.bmk"
    threads: MAX_THREADS
    conda: "homologs.yml"
    shell:
        """
        PATH="bin:$PATH"
        python src/pdc2/scripts/fasta_to_tree_pxclsq.py \
            {params.in_dir} \
            {threads} \
            aa \
            n \
        2> {log} 1>&2

        find {HOMOLOGS_R2} -name "OG*_RS_*.txt" -delete
        rename.ul .raxml.tre .tre {HOMOLOGS_R2}/OG*.raxml.tre
        """

rule homologs_round2_trim_tips:
    """Trim tips via the old method"""
    input:
        HOMOLOGS + "round2_fasta_to_tree.ok"
    output: touch(HOMOLOGS + "round2_trim_tips.ok")
    log: HOMOLOGS + "round2_trim_tips.log"
    benchmark: HOMOLOGS + "round2_trim_tips.bmk"
    params:
        relative_cutoff=params["homologs"]["trim_tips"]["relative_cutoff"],
        absolute_cutoff=params["homologs"]["trim_tips"]["absolute_cutoff"]
    conda: "homologs.yml"
    shell:
        """
        python src/pdc2/scripts/trim_tips.py \
            {HOMOLOGS_R2} \
            .tre \
            {params.relative_cutoff} \
            {params.absolute_cutoff} \
        2> {log} 1>&2
        """


rule homologs_round2_mask_tips_by_taxon_id:
    input: HOMOLOGS + "round2_trim_tips.ok"
    output: touch(HOMOLOGS + "round2_mask_tips_by_taxon_id.ok")
    params:
        in_dir = HOMOLOGS_R2,
        mask_tips = MASK_PARAPHYLETIC_TIPS
    threads: 1
    log: HOMOLOGS + "round2_mask_tips_by_taxon_id.log",
    benchmark: HOMOLOGS + "round2_mask_tips_by_taxon_id.bmk"
    conda: "homologs.yml"
    shell:
        "python src/pdc2/scripts/mask_tips_by_taxonID_transcripts.py "
            "{params.in_dir} "
            "{params.in_dir} "
            "{params.mask_tips} "
        "2> {log} 1>&2"


rule homologs_round2_cut_internal_long_branches:
    input: HOMOLOGS + "round2_mask_tips_by_taxon_id.ok"
    output: touch(HOMOLOGS + "round2_cut_internal_long_branches.ok")
    params:
        in_dir = HOMOLOGS_R2,
        internal_branch_cutoff = INTERNAL_BRANCH_CUTOFF,
        minimum_taxa = MINIMUM_TAXA
    threads: 1
    log: HOMOLOGS + "round2_cut_internal_long_branches.log"
    benchmark: HOMOLOGS + "homologs_round2_cut_internal_long_branches.bmk"
    conda: "homologs.yml"
    shell:
        "python src/pdc2/scripts/cut_long_internal_branches.py "
            "{params.in_dir} "
            ".mm "
            "{params.internal_branch_cutoff} "
            "{params.minimum_taxa} "
            "{params.in_dir} "
        "2> {log} 1>&2"


rule homologs_round2_write_fasta_files_from_trees:
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
        """
        python src/pdc2/scripts/write_fasta_files_from_trees.py \
            {input.fasta} \
            {params.indir} \
            .subtree \
            {params.indir} \
        2> {log} 1>&2

        find {HOMOLOGS_R2} -name "OG*rr.fa" -type f -exec \
            bash -c 'mv $1 ${{1%rr.fa}}.fa' _ {{}} \;
        """


rule homologs_round2:
    input:
        rules.homologs_round2_write_fasta_files_from_trees.output


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

rule homologs_rt_prepare:
    input: rules.homologs_round2.input
    output: touch(HOMOLOGS + "rt_prepare.ok")
    params:
        in_dir = HOMOLOGS_R2,
        out_dir = HOMOLOGS_RT
    log: HOMOLOGS + "rt_prepare.log"
    benchmark: HOMOLOGS + "rt_prepare.bmk"
    shell:
        """
        mkdir -p {params.out_dir}
        ln {params.in_dir}/*.subtree {params.out_dir}
        """


rule homologs_rt_prune_paralogs:
    input:
        tsv = HOMOLOGS + "in_out.tsv",
        ok = HOMOLOGS + "rt_prepare.ok"
    output:
        ok = touch(HOMOLOGS + "rt_prune_paralogs.ok")
    params:
        in_dir = HOMOLOGS_R2,
        tree_ending = ".subtree",
        out_dir = HOMOLOGS_RT,
        minimum_ingroup_taxa = MINIMUM_INGROUP_TAXA
    log: HOMOLOGS + "rt_prune_paralogs.log"
    benchmark: HOMOLOGS + "rt_prune_paralogs.bmk"
    conda: "homologs.yml"
    shell:
        """
        python src/pdc2/scripts/prune_paralogs_RT.py \
            {params.in_dir} \
            {params.tree_ending} \
            {params.out_dir} \
            {params.minimum_ingroup_taxa} \
            {input.tsv} \
        2> {log} 1>&2
        """

rule homologs_rt:
    input: HOMOLOGS + "rt_prune_paralogs.ok"


rule homologs_refine1_trees_to_fasta:
    input: 
        fasta = HOMOLOGS + "all.pep",
        ok = HOMOLOGS + "rt_prune_paralogs.ok"
    output: 
        touch(HOMOLOGS + "refine1_trees_to_fasta.ok"),
        directory(HOMOLOGS_REFINE1)
    params:
        indir = HOMOLOGS_RT
    log: HOMOLOGS + "refine1_trees_to_fasta.log"
    benchmark: HOMOLOGS + "refine1_trees_to_fasta.bmk"
    conda: "homologs.yml"
    shell:
        """
        python src/pdc2/scripts/write_fasta_files_from_trees.py \
            {input.fasta} \
            {params.indir} \
            .tre \
            {output[1]} \
        2> {log} 1>&2

        rename.ul rr.fa .fa {HOMOLOGS_REFINE1}/OG*rr.fa
        """


rule homologs_refine1:
    input: HOMOLOGS + "refine1_trees_to_fasta.ok"
    output: touch(HOMOLOGS + "refine1.ok")
    threads: MAX_THREADS
    log: HOMOLOGS + "refine1.log"
    benchmark: HOMOLOGS + "refine1.bmk"
    conda: "homologs.yml"
    shell:
        """
        python2 src/refine_alignments2.py \
            {HOMOLOGS_REFINE1} \
            .fa \
            {threads} \
        2> {log} 1>&2
        """


rule homologs_refine2_prepare:
    input: HOMOLOGS + "refine1.ok"
    output:
        touch(HOMOLOGS + "refine2_prepare.ok"),
        directory(HOMOLOGS_REFINE2)
    shell:
        """
        mkdir -p {HOMOLOGS_REFINE2}

        find {HOMOLOGS_REFINE1} -type f -name "*.maxalign.fa" \
        | parallel -j 1 basename {{}} .maxalign.fa \
        | parallel \
            ln --symbolic --relative \
                {HOMOLOGS_REFINE1}{{}}.maxalign.fa \
                {HOMOLOGS_REFINE2}{{}}.fa 
        """


rule homologs_refine2:
    input: HOMOLOGS + "refine2_prepare.ok"
    output: touch(HOMOLOGS + "refine2.ok")
    threads: MAX_THREADS
    log: HOMOLOGS + "refine2.log"
    benchmark: HOMOLOGS + "refine2.bmk"
    conda: "homologs.yml"
    shell:
        """
        python src/refine_alignments2.py \
            {HOMOLOGS_REFINE2} \
            .fa \
            {threads} \
        2> {log} 1>&2
        """


rule homologs:
    input: HOMOLOGS + "refine2.ok"