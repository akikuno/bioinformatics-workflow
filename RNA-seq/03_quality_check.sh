#!/bin/bash

###############################################################################
# Define threads
###############################################################################

[ -z "$threads" ] && threads=$(getconf _NPROCESSORS_ONLN 2>/dev/null | awk '{print int($0) - 1}')

###############################################################################
# FASTQC
###############################################################################

mkdir -p reports/fastqc

for file in data/fastq/*; do
    fastqc -t "$threads" "$file" \
        --nogroup -o reports/fastqc
done
