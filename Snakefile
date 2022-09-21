rule download_genomes:
    input:
        script: "code/build_kraken_database.sh"
    output:
