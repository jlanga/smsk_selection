rule download_uniprot_sprot:
    output:
        DOWNLOAD + "uniprot_sprot.dat.gz"
    params:
        url = features["swissprot"]
    log:
        DOWNLOAD + "uniprot_sprot.log"
    benchmark:
        DOWNLOAD + "uniprot_sprot.json"
    conda:
        "download.yml"
    shell:
        """
        wget \
            --continue \
            --output-document {output} \
            {params.url} \
        2> {log}
        """


rule download_nog_annotations:
    output:
        DOWNLOAD + "NOG.annotations.tsv.gz"
    params:
        url = features["NOG.annotations"]
    log:
        DOWNLOAD + "NOG.annotations.log"
    benchmark:
        DOWNLOAD + "NOG.annotations.json"
    conda:
        "download.yml"
    shell:
        """
        wget \
            --continue \
            --output-document {output} \
            {params.url} \
        2> {log}
        """


rule download_obo:
    output:
        DOWNLOAD + "go-basic.obo"
    params:
        url = features["obo"]
    log:
        DOWNLOAD + "obo.log"
    benchmark:
        DOWNLOAD + "obo.json"
    conda:
        "download.yml"
    shell:
        """
        wget \
            --continue \
            --output-document {output} \
            {params.url} \
        2> {log}
        """


rule download_pfama:
    output:
        DOWNLOAD + "Pfam-A.hmm.gz"
    params:
        url = features["Pfam-A"]
    log:
        DOWNLOAD + "pfama.log"
    benchmark:
        DOWNLOAD + "pfama.json"
    conda:
        "download.yml"
    shell:
        """
        wget \
            --continue \
            --output-document {output} \
            {params.url} \
        2> {log}
        """


rule download_busco_database:
    output:
        DOWNLOAD + "{database}_odb9.tar.gz"
    params:
        url = "http://busco.ezlab.org/v2/datasets/{database}_odb9.tar.gz"
    conda:
        "download.yml"
    log:
        DOWNLOAD + "busco_{database}_odb9.log"
    benchmark:
        DOWNLOAD + "busco_{database}_odb9.bmk"
    shell:
        """
        wget --continue {params.url} --output-document {output} 2> {log} 1>&2
        """
