rule download_wget_uniref90:
    output:
        fa_gz = protected(
            download + "uniref90.fa.gz"
        )
    threads:
        1
    params:
        url = config["databases"]["uniref90"]["url"]
    log:
        download + "wget_uniref90.log"
    benchmark:
        download + "wget_uniref90.json"
    shell:
        "wget "
            "--continue "
            "--output-document {output.fa_gz} "
            "{params.url} "
        "2> {log}"



rule download_wget_swissprot:
    output:
        fa_gz = protected(
            download + "swissprot.fa.gz"
        )
    threads:
        1
    params:
        url = config["databases"]["swissprot"]["url"]
    log:
        download + "wget_swissprot.log"
    benchmark:
        download + "wget_swissprot.json"
    shell:
        "wget "
            "--continue "
            "--output-document {output.fa_gz} "
            "{params.url} "
        "2> {log}"




rule download_wget_pfama:
    output:
        hmm_gz = protected(
            download + "pfama.hmm.gz"
        )
    threads:
        1
    params:
        url = config["databases"]["pfama"]["url"]
    log:
        download + "wget_pfama.log"
    benchmark:
        download + "wget_pfams.json"
    shell:
        "wget "
            "--continue "
            "--output-document {output.hmm_gz} "
            "{params.url} "
        "2> {log}"



