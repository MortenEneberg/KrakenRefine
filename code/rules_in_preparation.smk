


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
