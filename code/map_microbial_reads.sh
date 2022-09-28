#!/usr/bin/env bash

module purge

Bowtie2=$1
reads_taxID=$2
sam=$3
report=$4

module load $Bowtie2

bowtie2 -x data/index/kraken2_genomes \
	-U $reads_taxID \
	-f \
	--local \
	-D 30 -R 4 -N 0 -L 5 -i S,1,0.25 \
	-p 50 \
	-a \
	--mp 2,2 \
	--score-min C,31 \
	--met-file $report \
	-S $sam


module purge