configfile: "config.yaml"

rule all:
    input:
        expand("data/mapped_reads/{sample}/mapping_microbial_reads.done", sample = config["SAMPLE"])

rule build_bowtie2_database:
    output: 
        multiext("data/index/kraken2_genomes", ".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2"),
        touch("data/index/build_database.done")
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

rule load_database_to_shm:
    input:
        kraken2_database = config["kraken2db"]
    output:
        touch("data/kraken2_classification/database_loaded.done")
    shell:
        """
        mkdir /dev/shm/NCBI_GTDB_merge_kraken2db
        cp {input.kraken2_database}/*k2d /dev/shm/NCBI_GTDB_merge_kraken2db/
        """

rule classify_reads_kraken2:
    input: 
        reads = config["reads"]+"{sample}.fastq",
        loaded_database = "data/kraken2_classification/database_loaded.done"
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
        """

rule unload_database:
    input: 
        expand("data/kraken2_classification/{sample}/{sample}.report", sample = config["SAMPLE"])
    output:
        touch("data/kraken2_classification/database_unloaded.done")
    shell:
        """
        rm -r /dev/shm/NCBI_GTDB_merge_kraken2db
        """

rule extract_microbial_taxids:
    input:
        "data/kraken2_classification/{sample}/{sample}.report",
        config["Rlibpath"],
        "data/kraken2_classification/database_unloaded.done"
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
            shell('module load Biopython/1.78-foss-2020b-Python-3.8.6; python code/extract_kraken_reads.py -k {input.kraken_file} -s {input.reads} -o "data/microbial_reads/{wildcards.sample}/' + taxID + '.fastq" -t ' + taxID + ' -r {input.report_file} --include-children; module purge')


rule map_microbial_reads:
    input:
        "data/index/build_database.done",
        taxID_list="data/microbial_taxids/{sample}/genus_list.txt",
        done_statement_extract_reads = "data/microbial_reads/{sample}/extracting_microbial_reads.done"
    params:
        bowtie2 = config["Bowtie2"]
    output:
        touch("data/mapped_reads/{sample}/mapping_microbial_reads.done")
    run:
        import os
        genus_list = SplitGenusList(input.taxID_list)
        for taxID in genus_list:
            shell('code/map_microbial_reads.sh {params.bowtie2} "data/microbial_reads/{wildcards.sample}/' + taxID + '.fastq" "data/mapped_reads/{wildcards.sample}/' + taxID + '.sam" "data/mapped_reads/{wildcards.sample}/' + taxID + '.csv"') 

rule extract_read_lengths:
    input:
        reads = config["reads"]+"{sample}.fastq",
    output:
        read_lengths = "data/read_lengths/{sample}/read_lengths.tsv"
    params:
        seqkit = config["Seqkit"]
    shell:
        """    
        module load {params.seqkit}
        seqkit fx2tab -nl {input.reads} > {output.read_lengths}
        module purge
        """

rule evaluate_mappings:
    input:
        "data/mapped_reads/{sample}/mapping_microbial_reads.done",
         config["Rlibpath"],
        "data/kraken2_classification/{sample}.report",
        read_lengths = "data/read_lengths/{sample}/read_lengths.tsv"
    output:
        "data/tax_filtered_report/{sample}.report"
    script:
        "code/evaluate_mappings.R"

