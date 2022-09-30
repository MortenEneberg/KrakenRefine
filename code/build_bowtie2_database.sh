#!/usr/bin/env bash

module purge

Outdir=$3
Bowtie=$2
genomes=$1



module load $Bowtie

bowtie2-build -f $genomes $Outdir/kraken2_genomes --threads 20 -q

module purge
