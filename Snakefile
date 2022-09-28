configfile: "config.yaml"
#shell("""
#    mkdir /dev/shm/NCBI_GTDB_merge_kraken2db/
#    cp "/user_data/men/sepseq/databases/2022_05_25_new_kraken_database/kraken2_db_human_genome/*.k2d" /dev/shm/NCBI_GTDB_merge_kraken2db/
#""")

rule all:
    input:
        multiext("data/index/kraken2_genomes", ".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2"),
        expand("data/kraken2_classification/{sample}/{sample}.report", sample = config["SAMPLE"]),
        expand("data/microbial_taxids/{sample}/genus_list.txt", sample = config["SAMPLE"]),
        expand("data/microbial_reads/{sample}/extracting_microbial_reads.done", sample = config["SAMPLE"]),
        expand("data/mapped_reads/{sample}/mapping_microbial_reads.done", sample = config["SAMPLE"])


rule build_bowtie2_database:
    output: multiext("data/index/kraken2_genomes", ".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2")
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
        reads = config["reads"]+"{sample}.fastq"
    params:
        kraken2 = config["Kraken2"],
        kraken2_database = config["kraken2db"],
        confidence = config["Kraken2_confidence"],
        outdir = "data/kraken2_classification/{sample}/"
    output:
        report = "data/kraken2_classification/{sample}/{sample}.report",
        kraken = "data/kraken2_classification/{sample}/{sample}.kraken"
    shell:
        """
        code/classify_reads_kraken2.sh {params.kraken2} \
        {params.kraken2_database} \
        {params.confidence} \
        {params.outdir} \
        {input.reads} \
        {output.report} \
        {output.kraken}

        #rm -r /dev/shm/NCBI_GTDB_merge_kraken2db
        """


rule extract_microbial_taxids:
    input:
        "data/kraken2_classification/{sample}/{sample}.report",
        config["Rlibpath"]
    output:
        "data/microbial_taxids/{sample}/genus_list.txt"
    script:
        "code/extract_microbial_taxids.R"

def SplitGenusList(genus_list):
    genus_list = open(genus_list)
    contents = genus_list.read()
    contents_split = contents.splitlines()
    return contents_split

rule extract_microbial_reads:
    input:
        taxID_list="data/microbial_taxids/{sample}/genus_list.txt",
        kraken_file="data/kraken2_classification/{sample}/{sample}.kraken",
        reads=config["reads"] + "{sample}.fastq",
        report_file="data/kraken2_classification/{sample}/{sample}.report"
    output:
        touch("data/microbial_reads/{sample}/extracting_microbial_reads.done")
    run:
        import os
        genus_list = SplitGenusList(input.taxID_list)
        for taxID in genus_list:
            shell('python code/extract_kraken_reads.py -k {input.kraken_file} -s {input.reads} -o "data/microbial_reads/{wildcards.sample}/' + taxID + '_.fastq" -t ' + taxID + ' -r {input.report_file} --include-children')


rule map_microbial_reads:
    input:
        reads_taxID = "data/microbial_reads/{sample}/{taxID}.fastq"#,
        #multiext("data/index/kraken2_genomes.1.bt2", ".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2")
    params:
        bowtie2 = config["Bowtie2"]
    output:
        #sam = "data/mapped_reads/{sample}/{taxID}.sam",
        #report = "data/mapped_reads/{sample}/{taxID}.csv",
        touch("data/mapped_reads/{sample}/mapping_microbial_reads.done")
    shell:
        """
        code/map_microbial_reads.sh \
        {params.bowtie2} \
        {input.reads_taxID} \
        {output.sam} \
        {output.report}
        """