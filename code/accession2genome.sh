#!/usr/bin/env bash

module purge

output=$1
genomes=$2

find $genomes -type f -name "*.fna" -exec grep -rH ">" {} \; | sed -r 's/[:>]+/\t/g' > $output