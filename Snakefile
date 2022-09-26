configfile: "config.yaml"

rule all:
    input:
        "data/index/kraken2_genomes.1.bt2",
        expand("data/kraken2_classification/{sample}.report", sample = config["SAMPLE"])


rule build_bowtie2_database:
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

rule load_db_to_shm:
    output:
        "dev/shm/NCBI_GTDB_merge_kraken2db/hash.k2d"
    params:
        kraken2db = config["kraken2db"]
    shell:
        """
        mkdir /dev/shm/NCBI_GTDB_merge_kraken2db/
        cp {params.kraken2db}*k2d /dev/shm/NCBI_GTDB_merge_kraken2db/
        """

rule classify_reads_kraken2:
    input: 
        "dev/shm/NCBI_GTDB_merge_kraken2db/hash.k2d",
        reads = config["reads"]+"{sample}.fastq"
    params:
        kraken2 = config["Kraken2"],
        kraken2_database = config["kraken2db"],
        confidence = config["Kraken2_confidence"],
        outdir = "data/kraken2_classification"
    output:
        report = "data/kraken2_classification/{sample}.report",
        kraken = "data/kraken2_classification/{sample}.kraken"
    shell:
        """
        code/classify_reads_kraken2.sh \ 
        {params.kraken2} \
        {params.kraken2_database} \
        {params.confidence} \
        {params.outdir} \
        {input.reads} \
        {output.report} \
        {output.kraken}

        rm -r /dev/shm/NCBI_GTDB_merge_kraken2db
        """