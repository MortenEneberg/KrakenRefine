#! /bin/bash

## Build a database from the genomes also used in the Kraken2 database to which all non-human reads will be aligned to.

#module load Bowtie2/2.4.2-foss-2020b
#bowtie2-build -f $sample_dir/$sample_name.fasta $sample_name --threads 20 

echo "Hi"

#$snakemake_config[reads]