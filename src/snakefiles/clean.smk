rule clean_raw:
    shell:
        "rm -ri "
      		"{raw_dir}"


rule clean:
    shell:
        "rm -rf results"
