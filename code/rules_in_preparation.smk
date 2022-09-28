
rule evaluate_mappings:
    input:
        "data/mapped_reads/{sample}/{taxID}.sam",
        "data/kraken2_classification/{sample}.report"
    output:
        "data/tax_filtered_report/{sample}.report"
    script:
        "code/evaluate_mappings.R"

