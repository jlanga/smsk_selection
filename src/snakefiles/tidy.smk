rule tidy_pep:
    input: CDHIT + "{species}.pep"
    output: TIDY + "{species}.pep"
    params: "{species}_"
    conda: "tidy.yml"
    shell:
        """
        < {input} \
        bash src/tag_fasta.sh {params} \
        | bash src/remove_stops_pep.sh \
        > {output}
        """

rule tidy_cds:
    input: CDHIT + "{species}.cds"
    output: TIDY + "{species}.cds"
    params: "{species}_"
    conda: "tidy.yml"
    shell:
        """
        < {input} \
        bash src/tag_fasta.sh {params} \
        | bash src/remove_stops_cds.sh \
        > {output}
        """