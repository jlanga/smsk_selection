rule filterlen_pep:
    input: TAG + "{species}.pep"
    output: FILTERLEN + "{species}.pep"
    conda: "filterlen.yml"
    shell: "python src/filter_longest_protein_per_gene.py {output} {input}"


rule filterlen_cds:
    input:
        tag_cds = TAG + "{species}.cds",
        pep_fai = FILTERLEN + "{species}.pep.fai"
    output: FILTERLEN + "{species}.cds"
    conda: "filterlen.yml"
    shell:
        """cut -f 1 {input.pep_fai} \
        | cut -f 1 -d . \
        | xargs seqtk seq {input.tag.cds} \
        > {output}
        """

rule filterlen:
    input:
        expand(
            FILTERLEN + "{species}.{extension}",
            species=SPECIES,
            extension="pep cds".split()
        )
