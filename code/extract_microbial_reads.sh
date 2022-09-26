#!/usr/bin/env bash

module purge

Kraken2=$1
kraken2db=$2
confidence=$3
Outdir=$4
reads=$5
output_report=$6
output_kraken=$7

mkdir $4




rule extract_microbial_taxids:
    input:
        "data/kraken2_classification/{sample}.report",
        config["Rlibpath"]
    output:
        "data/microbial_taxids/{sample}/genus_{taxID}.txt"
    script:
        "code/extract_microbial_taxids.R"

rule extract_microbial_reads:
    input:
        taxid_lists = "data/{sample}_genus_taxids.txt",
        kraken_files = "data/kraken2_classification/{sample}.kraken",
        reads = config["reads"]+"{sample}.fastq"
    params:

    output:
        "data/microbial_reads/{sample}/{taxID}.fastq"
    run:
        f = open({input.taxid_lists},'r')


        code/extract_kraken_reads.py


rule map_microbial_reads:
    input:
    output:
        "data/mapped_reads/{sample}/{taxID}.sam"
    shell:
        """
        code/map_microbial_reads.sh \
        {params.output}
        """

rule evaluate_mappings:
    input:
        "data/mapped_reads/{sample}/{taxID}.sam"
    output:
        "data/map_evaluation/{sample}/taxID_evaluated.txt" ##Could be a TAB separated file with a taxID and a present or not-present flag to it
    script:
        "code/evaluate_mappings.R"


