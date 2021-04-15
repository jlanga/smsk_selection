def get_species_pep(wildcards):
    return samples.loc[wildcards.species]["pep"]


def get_species_cds(wildcards):
    return samples.loc[wildcards.species]["cds"]


rule raw_link_pep:
    output: RAW + "{species}.pep"
    params: get_species_pep
    shell: "ln --symbolic $(readlink --canonicalize {params}) {output}"


rule raw_link_cds:
    output: RAW + "{species}.cds"
    params: get_species_cds
    shell: "ln --symbolic $(readlink --canonicalize {params}) {output}"
