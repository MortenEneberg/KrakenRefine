


rule map_microbial_reads:
    input:
        reads_taxID = "data/microbial_reads/{sample}/{taxID}.fastq",
        multiext("data/index/kraken2_genomes", ".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2")
    params:
        bowtie2 = config["Bowtie2"]
    output:
        sam = "data/mapped_reads/{sample}/{taxID}.sam",
        report = "data/mapped_reads/{sample}/{taxID}.csv"
    shell:
        """
        code/map_microbial_reads.sh \
        {params.bowtie2} \
        {input.reads_taxID} \
        {output.sam} \
        {output.report}
        """

rule evaluate_mappings:
    input:
        "data/mapped_reads/{sample}/{taxID}.sam",
        "data/kraken2_classification/{sample}.report"
    output:
        "data/tax_filtered_report/{sample}.report"
    script:
        "code/evaluate_mappings.R"

