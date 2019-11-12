rule orthogroups_join_pep:
    """
    Put together all .pep files into all.pep
    """
    input:
        pep = expand(
            tag + "{species}.pep",
            species = species
        )
    output:
        pep = orthogroups + "all.pep"
    threads:
        1
    log:
        orthogroups + "join_pep.log"
    benchmark:
        orthogroups + "join_pep.json"
    shell:
        "cat {input.pep} > {output.pep} 2> {log}"



rule orthogroups_join_cds:
    """
    Put together all .cds files into all.cds
    """
    input:
        cds = expand(
            tag + "{species}.cds",
            species = species
        )
    output:
        cds = orthogroups + "all.cds"
    threads:
        1
    log:
        orthogroups + "join_cds.log"
    benchmark:
        orthogroups + "join_cds.json"
    shell:
        "cat {input.cds} > {output.cds} 2> {log}"



rule orthogroups_extract_pep:
    """
    Extract into a fasta file the sequences of each orthogroup
    """
    input:
        pep = orthogroups + "all.pep",
        groups_txt = orthofinder + "OrthologousGroups.txt"
    output:
        dynamic(
            orthogroups_pep + "{orthogroup_id}.fasta"
        )
    params:
        folder = orthogroups_pep,
        minimum = config["params"]["orthogroups"]["extract_orthogroups"]["minimum"]
    threads:
        1
    log:
        orthogroups + "extract_pep.log"
    benchmark:
        orthogroups + "extract_pep.json"
    shell:
        "python3 bin/extract_single_copy_orthologs.py "
            "--out_dir {params.folder} "
            "--infasta {input.pep} "
            "--groups {input.groups_txt} "
            "--species {species} "
        "2> {log} 1>&2"



rule orthogroups_remove_stops_pep:
    """
    Remove stops (*) at the end of each sequence in a fasta file (pep)
    """
    input:
        fasta = orthogroups_pep + "{orthogroup_id}.fasta"
    output:
        fasta = orthogroups_pep_nonstop + "{orthogroup_id}.fasta"
    threads:
        1
    log:
        orthogroups_pep + "log/remove_stops_{orthogroup_id}.log"
    benchmark:
        orthogroups_pep + "log/remove_stops_{orthogroup_id}.json"
    run:
        with open(input.fasta, "r") as f_in, open(output.fasta, "w") as f_out:
            for line in f_in:
                f_out.write(line.replace("*\n", "\n"))





rule orthogroups_align_pep:
    input:
        fasta = orthogroups_pep_nonstop + "{orthogroup_id}.fasta"
    output:
        fasta = orthogroups_pep_msa + "{orthogroup_id}.fasta"
    threads:
        1
    log:
        orthogroups + "log/align_{orthogroup_id}.log"
    benchmark:
        orthogroups + "log/align_{orthogroup_id}.json"
    shell:
        "mafft "
            "--auto "
            "--thread 1 "
            "--quiet "
            "{input.fasta} "
        "> {output.fasta} "
        "2> {log}"



rule orthogroups_trim_pep:
    input:
        fasta = orthogroups_pep_msa + "{orthogroup_id}.fasta"
    output:
        fasta = orthogroups_pep_trim + "{orthogroup_id}.fasta"
    threads:
        1
    log:
        orthogroups + "log/trim_pep_{orthogroup_id}.log"
    benchmark:
        orthogroups + "log/trim_pep_{orthogroup_id}.json"
    shell:
        "trimal "
            "-in {input.fasta} "
            "-out {output.fasta} "
            "-automated1 "
        "2> {log} 1>&2"



rule orthogroups_extract_cds:
    """
    Extract into a fasta file the sequences of each orthogroup
    """
    input:
        cds = orthogroups + "all.cds", 
        groups_txt = orthofinder + "OrthologousGroups.txt"
    output:
        dynamic(
            orthogroups_cds + "{orthogroup_id}.fasta"
        )
    params:
        folder = orthogroups_cds,
        minimum = config["params"]["orthogroups"]["extract_orthogroups"]["minimum"]
    threads:
        1
    log:
        orthogroups + "extract_cds.log"
    benchmark:
        orthogroups + "extract_cds.json"
    shell:
        "python3 bin/extract_single_copy_orthologs.py "
            "--out_dir {params.folder} "
            "--infasta {input.cds} "
            "--groups {input.groups_txt} "
            "--species {species} "
        "2> {log} 1>&2"



rule orthogroups_pal2nal_orthogroup:
    """
    Convert protein alignment into phy dna alignment
    """
    input:
        pep_aln = orthogroups_pep_msa + "{orthogroup_id}.fasta",
        cds = orthogroups_cds + "{orthogroup_id}.fasta"
    output:
        phy = orthogroups_cds_msa + "{orthogroup_id}.phy"
    threads:
        1
    log:
        orthogroups + "log/pal2nal_{orthogroup_id}.log"
    benchmark:
        orthogroups + "log/pal2nal_{orthogroup_id}.json"
    shell:
        "( perl src/pal2nal.v14/pal2nal.pl "
            "{input.pep_aln} "
            "{input.cds} "
            "-output fasta | "
        "python3 bin/fasta_to_phy.py "
        "> {output} ) "
        "2> {log}"



rule orthogroups_phyml_orthogroup:
    """
    Get the ML tree with Phyml
    """
    input:
        phy = orthogroups_cds_msa + "{orthogroup_id}.phy"
    output:
        nwk = orthogroups + "phyml/{orthogroup_id}.nwk",
        stats = orthogroups + "phyml/log/{orthogroup_id}.stats"
    threads:
        1
    log:
        orthogroups + "phyml/log/{orthogroup_id}.log"
    benchmark:
        orthogroups + "phyml/log/{orthogroup_id}.json"
    shell:
        "( phyml "
            "-i {input.phy} "
            "-d aa ; "
        "mv "
            "{input.phy}_phyml_tree "
            "{output.nwk} ; "
        "mv "
            "{input.phy}_phyml_stats "
            "{output.stats} ) "
        "2> {log} 1>&2"



rule orthogroups_fastcodeml_orthogroup:
    """
    Get patterns of positive selection with fastcodeml
    """
    input:
        nwk = orthogroups + "phyml/{orthogroup_id}.nwk",
        phy = orthogroups_cds_msa + "{orthogroup_id}.phy"
    output:
        fastcodeml = orthogroups + "fastcodeml/{orthogroup_id}.fastcodeml"
    threads:
        1
    params:
        maximizer = config["params"]["orthogroups"]["fastcodeml"]["maximizer"],
    log:
        orthogroups + "fastcodeml/log/{orthogroup_id}.log"
    benchmark:
        orthogroups + "fastcodeml/log/{orthogroup_id}.json"
    shell:
        "./bin/fast "
            "--maximizer {params.maximizer} "
            "--branch-from-file "
            "{input.nwk} "
            "{input.phy} "
        "> {output.fastcodeml} "
        "2> {log}"

