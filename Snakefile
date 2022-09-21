configfile: "config.yaml"


rule build_bowtie_database:
    output: "data/index/kraken2_genomes.1.bt2"
    threads: 30
    params: 
        genomes = config["kraken2_genomes"],
        bowtie2 = config["Bowtie2"],
        outdir = "data/index"
    shell:
        """
        code/build_bowtie2_database.sh \
        {params.genomes} \
        {params.bowtie2} \
        {params.outdir}
        """


rule classify_reads_kraken2:
    input: 
        reads = config["reads"]
    params:
        kraken2 = config["Kraken2"],
        kraken2_database = config["kraken2db"],
        confidence = config["Kraken2_confidence"],
        outdir = "data/kraken2_classification"
    output:
        
    shell:
        """
        code/classify_reads_kraken2.sh \ 
        {params.kraken2} \
        {params.kraken2_database} \
        {params.confidence} \
        {params.outdir}
        """