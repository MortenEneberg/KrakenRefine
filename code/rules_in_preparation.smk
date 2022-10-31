
rule evaluate_mappings:
    input:
        touch("data/mapped_reads/{sample}/mapping_microbial_reads.done"),
         config["Rlibpath"],
        "data/kraken2_classification/{sample}.report",
        read_lengths = "data/read_lengths/{sample}/read_lengths.tsv"
    output:
        "data/tax_filtered_report/{sample}.report"
    script:
        "code/evaluate_mappings.R"


gunzip -r "/user_data/men/sepseq/databases/2022_10_18_kraken2_EUPATH_database/genomes/"

grep -r ">" *.fna | sed -r 's/[:>]+/\t/g' > /user_data/men/sepseq/clinical_studies/init_tests/kraken2_filter/data/mapped_reads/accession2genome.tsv
