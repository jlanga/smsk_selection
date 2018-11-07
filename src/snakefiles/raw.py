def get_species_fasta(wildcards):
    return samples.loc[wildcards.species]["fasta"]


def get_species_g2t(wildcards):
    return samples.loc[wildcards.species]["g2t"]


def get_species_url(wildcards):
    return samples.loc[wildcards.species]["url"]


rule raw_link_fa:
    output: RAW + "{species}.fasta"
    params: get_species_fasta
    shell: "ln --symbolic $(readlink --canonicalize {params}) {output}"


rule raw_link_g2t:
    output: RAW + "{species}.g2t.tsv"
    params: get_species_g2t
    shell: "ln --symbolic $(readlink --canonicalize {params}) {output}"


rule raw_download_species:
    output: protected( "data/transcriptomes/" + "{species}.fasta")
    params: get_species_url
    log: RAW + "download_{species}.log"
    shell:
        """
        (wget --output-document - {params.url} \
        | pigz --decompress --stdout \
        > {output}) 2> {log}
        """
