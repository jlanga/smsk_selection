rule busco_run:
    input:
        fasta = RAW + "{species}.fasta",
        db_folder = DB + "{database}_odb9"
    output:
        directory(BUSCO + "{species}_{database}")
    log:
        BUSCO + "{species}_{database}.log"
    benchmark:
        BUSCO + "{species}_{database}.bmk"
    threads:
        4
    params:
        output_tag = "{species}_{database}"
    conda:
        "busco.yml"
    shell:
        """
        run_busco \
            --in {input.fasta} \
            --out {params.output_tag} \
            --lineage_path {input.db_folder} \
            --mode tran \
            --cpu {threads} \
        2> {log} 1>&2

        mv run_{params.output_tag} {output}
        """


rule busco:
    input:
        expand(
            BUSCO + "{species}_{database}",
            species=SPECIES,
            database=features["busco_species"]
        )
