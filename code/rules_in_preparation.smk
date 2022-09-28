


def SplitGenusList(genus_list):
    genus_list = open(genus_list)
    contents = genus_list.read()
    contents_split = contents.splitlines()
    return contents_split

checkpoint split_taxids:
    input:
        taxID_fastq=expand(result("data/microbial_reads/{sample}/{taxID}.fastq"), taxID=SplitGenusList("data/microbial_taxids/{sample}/genus_list.txt"))
    output:
        touch(result("taxID_introduction.done"))

taxIDs = None
def find_bins_with_16s():
    output_list=[]
    for name in SAMPLES:
        path = "cmscan_subsetted_16s/"+name+".tblout"
        num_extracted = helper_functions.count_number_of_lines_in_file(path)-1
        if num_extracted is not None and int(num_extracted) > 0:
            output_list.append(name)
    global taxIDs
    taxIDs = output_list


def get_bins_with_extracted_16s_to_combine(wildcards):
    checkpoint_output = checkpoints.collect_bins_with_16s.get()
    if taxIDs is None:
        find_bins_with_16s()
    names = taxIDs
    return expand("extracted_16s_with_mag_ids/{sample}.fa", sample=names)



rule extract_microbial_reads:
    input:
        taxID_list="data/microbial_taxids/{sample}/genus_list.txt",
        kraken_file="data/kraken2_classification/{sample}/{sample}.kraken",
        reads=config["reads"] + "{sample}.fastq"
    output:
        touch("data/microbial_reads/{sample}/extracting_microbial_reads.done")
    run:
        import os
        genus_list = SplitGenusList("data/microbial_taxids/{wildcards.sample}/genus_list.txt")
        for taxID in genus_list:
            os.system('python code/extract_kraken_reads.py --kraken {input.kraken_file} -s {input.reads} '
                      '-o "data/microbial_reads/{wildcards.sample}/' + taxID + '_.fastq" '
                      '-t ' + taxID + ' --include-children')

rule extract_microbial_reads:
    input:
        taxID_list = "data/microbial_taxids/{sample}/genus_list.txt",
        kraken_file = "data/kraken2_classification/{sample}/{sample}.kraken",
        reads = config["reads"]+"{sample}.fastq"
        genus_list = SplitGenusList(input.taxID_list)
    output:
        out_reads = expand("data/microbial_reads/{sample}/{taxID}.fastq", taxID = SplitGenusList(input.taxID_list))
    shell:
        """
        for taxID in genus_list
        do
        python code/extract_kraken_reads.py --kraken {input.kraken_file} -s {input.reads} -o "data/microbial_reads/{sample}/$taxID\_.fastq" -t $taxID --include-children
        done
        """


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

