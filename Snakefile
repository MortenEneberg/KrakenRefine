configfile: "config.yaml"
shell("""
    mkdir /dev/shm/NCBI_GTDB_merge_kraken2db/
    cp config["kraken2db"]*k2d /dev/shm/NCBI_GTDB_merge_kraken2db/
""")

rule all:
    input:
        multiext("data/index/kraken2_genomes", ".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2"),
        expand("data/kraken2_classification/{sample}/{sample}.report", sample = config["SAMPLE"]),
        expand("data/microbial_taxids/{sample}/genus_list.txt", sample = config["SAMPLE"])#,
   #     expand("data/microbial_reads/{sample}/{taxID}.fastq", taxID = SplitGenusList("data/microbial_taxids/{sample}/genus_list.txt"))


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

checkpoint split_taxids:
    input:
        taxID_fastq=expand(result("data/microbial_reads/{sample}/{taxID}.fastq"), taxID=SplitGenusList("data/microbial_taxids/{sample}/genus_list.txt"))
    output:
        touch(result("taxID_introduction.done"))

taxIDs = None
def find_bins_with_16s():
    output_list=[]
    for name in SAMPLES:
        path = result("cmscan_subsetted_16s/"+name+".tblout")
        num_extracted = helper_functions.count_number_of_lines_in_file(path)-1
        if num_extracted is not None and int(num_extracted) > 0:
            output_list.append(name)
    global BINS_WITH_16S
    BINS_WITH_16S = output_list
def get_bins_with_extracted_16s_to_combine(wildcards):
    checkpoint_output = checkpoints.collect_bins_with_16s.get()
    if BINS_WITH_16S is None:
        find_bins_with_16s()
    names = BINS_WITH_16S
    return expand(result("extracted_16s_with_mag_ids/{sample}.fa"), sample=names)



rule extract_microbial_reads:
    input:
        taxID_list = "data/microbial_taxids/{sample}/genus_list.txt",
        kraken_file = "data/kraken2_classification/{sample}/{sample}.kraken",
        reads = config["reads"]+"{sample}.fastq"
    output:
        out_reads = expand("data/microbial_reads/{sample}/{taxID}.fastq", taxID = SplitGenusList(input.taxID_list))
    run:
        for taxID in SplitGenusList({input.taxID_list}):
            shell("python code/extract_kraken_reads.py --kraken {input.kraken_file} -s {input.reads} -o {output.out_reads} -t taxID --include-children")

