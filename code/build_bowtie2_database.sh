#!/usr/bin/env bash

module purge

Outdir_index=$3
Outdir_genomes=$4
Bowtie=$2
genomes=$1


find $1 -name "*.fna.gz" -type f -exec cat {} + > $Outdir_genomes/cat_genomes.fna.gz
gunzip $Outdir_genomes/cat_genomes.fna.gz

module load $Bowtie

bowtie2-build -f $Outdir_genomes/cat_genomes.fna $Outdir_index/kraken2_genomes --threads 50 --large-index -q

module purge
