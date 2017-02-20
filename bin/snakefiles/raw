rule raw_make_link_species:
    output:
        fasta= raw + "{species}.fasta"
    threads:
        1
    log:
        raw + "make_link_{species}.log"
    benchmark:
        raw + "make_link_{species}.json"
    params:
        original_file= lambda wildcards: config["species"][wildcards.species]["fasta"]
    shell:
        "ln -s "
            "$(readlink -f {params.original_file}) "
            "{output.fasta} "
        "2> {log} "
        "1>&2"



rule raw_dowload_species:
    output:
        fasta = protected(
            "data/transcriptomes/" + "{species}.fasta"
        )
    params:
        url = lambda wildcards: config["species"][wildcards.species]["url"]
    threads:
        1
    log:
        raw + "download_{species}.log"
    benchmark:
        raw + "download_{species}.json"
    shell:
        "( wget "
            "--output-document - "
            "{params.url} | "
        "pigz "
            "--decompress "
            "--stdout "
        "> {output.fasta} ) "
        "2> {log}"
