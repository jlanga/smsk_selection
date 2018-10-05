rule filterlen_pep:
    input: tag + "{species}.pep"
    output: filterlen + "{species}.pep"
    conda: "filterlen.yml"
    shell: "python src/filter_longest_protein_per_gene.py {output} {input}"


rule filterlen_cds:
    input:
        tag_cds = tag + "{species}.cds",
        pep_fai = filterlen + "{species}.pep.fai"
    output: filterlen + "{species}.cds"
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
            filterlen + "{species}.{extension}",
            species=SPECIES,
            extension="pep cds".split()
        )
