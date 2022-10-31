#!/usr/bin/env bash

module purge

Bowtie2=$1
reads_taxID=$2
sam_no_header=$3
report_no_header=$4
sam_header=$5
report_header=$6

module load $Bowtie2

echo $reads_taxID

bowtie2 -x data/index/kraken2_genomes \
	-U $reads_taxID \
	-f \
	-p 20 \
	--local \
	--sensitive-local \
	-k 2 \
	--mp 2,2 \
	--score-min G,20,8 \
	--met-file $report_no_header \
	--no-hd \
	--mm \
	-S $sam_no_header

bowtie2 -x data/index/kraken2_genomes \
	-U $reads_taxID \
	-f \
	-p 20 \
	--local \
	--sensitive-local \
	-k 2 \
	--mp 2,2 \
	--score-min G,20,8 \
	--met-file $report_header \
	--mm \
	-S $sam_header

module purge