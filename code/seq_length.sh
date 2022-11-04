#!/usr/bin/env bash

module purge

output=$1
Samtools=$2
genomes=$3

module load $Samtools

samtools faidx $genomes 

cut -f1-2 data/genomes/cat_genomes.fna.fai > $output 

module purge
