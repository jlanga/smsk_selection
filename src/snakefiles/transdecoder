rule transdecoder_longorfs_species:
    """
    Predict ORFs by length
    """
    input:
        fasta = raw + "{species}.fasta"
    output:
        base_freqs = temp(
            "{species}.fasta.transdecoder_dir/base_freqs.dat"
        ),
        base_freqs_ok = temp(
            "{species}.fasta.transdecoder_dir/base_freqs.dat.ok"
        ),
        cds = temp(
            "{species}.fasta.transdecoder_dir/longest_orfs.cds"
        ),
        gff3 = temp(
            "{species}.fasta.transdecoder_dir/longest_orfs.gff3"
        ),
        pep = temp(
            "{species}.fasta.transdecoder_dir/longest_orfs.pep"
        ),
        dir = temp(
            "{species}.fasta.transdecoder_dir/"
        )
    threads:
        1
    log:
        transdecoder + "{species}/longorfs.log"
    benchmark:
        transdecoder + "{species}/longorfs.json"
    shell:
        "TransDecoder.LongOrfs "
            "-t {input.fasta} "
        "2> {log} 1>&2"



rule transdecoder_index_longest_orfs_species:
    """
    Create samtools index for longest_orfs.pep
    """
    input:
        pep = "{species}.fasta.transdecoder_dir/longest_orfs.pep"
    output:
        fai = "{species}.fasta.transdecoder_dir/longest_orfs.pep.fai"
    log:
        transdecoder + "{species}/index_longest_orfs.log"
    benchmark:
        transdecoder + "{species}/index_longest_orfs.json"
    shell:
        "samtools faidx {input.pep} 2> {log}"



rule transdecoder_split_longest_orfs_species:
    """
    Split the headers from transdecoder_longest_orfs into multiple files
    """
    input:
        fai = "{species}.fasta.transdecoder_dir/longest_orfs.pep.fai"
    output:
        expand(
            transdecoder + "{{species}}/chunks/longest_orfs_{chunk_id}.tsv",
            chunk_id = ['{0:05d}'.format(x) for x in range(0, config["params"]["transdecoder"]["number_of_chunks"])]
        )
    params:
        number_of_chunks = config["params"]["transdecoder"]["number_of_chunks"],
        species = "{species}"
    log:
        transdecoder + "{species}/chunks/split_longest_orfs.log"
    benchmark:
        transdecoder + "{species}/chunks/split_longest_orfs.json"
    shell:
        "split "
            "--number l/{params.number_of_chunks} "
            "--numeric-suffixes "
            "--suffix-length 5 "
            "--additional-suffix .tsv "
            "{input.fai} "
            "{transdecoder}/{params.species}/chunks/longest_orfs_ "
        "2> {log}"



rule transdecoder_hmmscan_species_chunk:
    """
    hmmscan over one chunk
    """
    input:
        pep = "{species}.fasta.transdecoder_dir/longest_orfs.pep",
        fai = "{species}.fasta.transdecoder_dir/longest_orfs.pep.fai",
        chunk = transdecoder + "{species}/chunks/longest_orfs_{chunk_id}.tsv",
        hmm = db + "pfama.hmm"
    output:
        tsv = transdecoder + "{species}/hmmscan/longest_orfs_{chunk_id}.tsv"
    log:
        transdecoder + "{species}/hmmscan/longest_orfs_{chunk_id}.log"
    benchmark:
        transdecoder + "{species}/hmmscan/longest_orfs_{chunk_id}.json"
    shell:
        "cut -f 1 {input.chunk} "
        "| xargs samtools faidx {input.pep} "
        "| hmmscan "
            "--domtblout {output.tsv} "
            "{input.hmm} "
            "- "
        "2> {log} 1>&2"



rule transdecoder_hmmscan_merge_species:
    """
    Merge hmmscan results into one file
    """
    input:
        expand(
            transdecoder + "{{species}}/hmmscan/longest_orfs_{chunk_id}.tsv",
            chunk_id = ['{0:05d}'.format(x) for x in range(0, config["params"]["transdecoder"]["number_of_chunks"])]
        )
    output:
        tsv = transdecoder + "{species}/hmmscan.tsv"
    log:
        transdecoder + "{species}/hmmscan_merge.log"
    benchmark:
        transdecoder + "{species}/hmmscan_merge.json"
    shell:
        "cat {input} > {output} 2> {log}"



rule transdecoder_blastp_species_chunk:
    """
    Run blastp of each chunk
    """
    input:
        pep = "{species}.fasta.transdecoder_dir/longest_orfs.pep",
        fai = "{species}.fasta.transdecoder_dir/longest_orfs.pep.fai",
        chunk = transdecoder + "{species}/chunks/longest_orfs_{chunk_id}.tsv",
        db = db + config["params"]["transdecoder"]["blastp"]["database"]
    output:
        tsv = transdecoder + "{species}/blastp/longest_orfs_{chunk_id}.tsv"
    log:
        transdecoder + "{species}/blastp/longest_orfs_{chunk_id}.log"
    benchmark:
        transdecoder + "{species}/blastp/longest_orfs_{chunk_id}.json"
    shell:
        "cut -f 1 {input.chunk} "
        "| xargs samtools faidx {input.pep} "
        "| blastp "
            "-db {input.db} "
            "-max_target_seqs 1 "
            "-outfmt 6 "
            "-evalue 1e-5 "
            "-out {output.tsv} " 
        "2> {log} 1>&2"



rule transdecoder_blastp_species_merge:
    """
    Merge results from the different blastps
    """
    input:
        expand(
            transdecoder + "{{species}}/blastp/longest_orfs_{chunk_id}.tsv",
            chunk_id = ['{0:05d}'.format(x) for x in range(0, config["params"]["transdecoder"]["number_of_chunks"])]
        )
    output:
        tsv = transdecoder + "{species}/blastp.tsv"
    log:
        transdecoder + "{species}/blastp_merge.log"
    benchmark:
        transdecoder + "{species}/blastp_merge.json"
    shell:
        "cat {input} > {output} 2> {log}"


rule transdecoder_predict_species:
    input:
        fasta = raw + "{species}.fasta",
        pfam_tsv = transdecoder + "{species}/hmmscan.tsv",
        blastp_tsv = transdecoder + "{species}/blastp.tsv",
        folder = "{species}.fasta.transdecoder_dir/",
        base_freqs = "{species}.fasta.transdecoder_dir/base_freqs.dat",
        base_freqs_ok = "{species}.fasta.transdecoder_dir/base_freqs.dat.ok",
        cds = "{species}.fasta.transdecoder_dir/longest_orfs.cds",
        gff3 = "{species}.fasta.transdecoder_dir/longest_orfs.gff3",
        pep = "{species}.fasta.transdecoder_dir/longest_orfs.pep"
    output:
        bed  = protected(transdecoder + "{species}.bed"),
        cds  = protected(transdecoder + "{species}.cds"),
        gff3 = protected(transdecoder + "{species}.gff3"),
        pep  = protected(transdecoder + "{species}.pep"),
    params:
        dir = "{species}.fasta.transdecoder_dir/",
        bed  = "{species}.fasta.transdecoder.bed",
        cds  = "{species}.fasta.transdecoder.cds",
        gff3 = "{species}.fasta.transdecoder.gff3",
        pep  = "{species}.fasta.transdecoder.pep",
    threads:
        24
    log:
        transdecoder + "predict_{species}.log"
    benchmark:
        transdecoder + "predict_{species}.json"
    shell:
        "TransDecoder.Predict "
            "-t {input.fasta} "
            "--retain_pfam_hits {input.pfam_tsv} "
            "--retain_blastp_hits {input.pfam_tsv} "
            "--cpu {threads} "
        "2> {log} 1>&2 && "
        "mv {params.bed} {output.bed} 2>> {log} 1>&2 && "
        "mv {params.cds} {output.cds} 2>> {log} 1>&2 && "
        "mv {params.gff3} {output.gff3} 2>> {log} 1>&2 && "
        "mv {params.pep} {output.pep} 2>> {log} 1>&2 && "
        "rm -rf {params.dir} 2>> {log} 1>&2"
