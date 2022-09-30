#!/usr/bin/env bash

module purge

Kraken2=$1
kraken2db=$2
confidence=$3
Outdir=$4
reads=$5
output_report=$6
output_kraken=$7


module load $Kraken2


kraken2 \
--db $kraken2db \
--confidence $confidence \
--memory-mapping \
--report $output_report \
--use-names \
--threads 20 --report-minimizer-data \
--output $output_kraken \
$reads


