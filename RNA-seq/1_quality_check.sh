#!/bin/bash

###############################################################################
# Define threads
###############################################################################

[ -z "$threads" ] && threads=$(getconf _NPROCESSORS_ONLN 2>/dev/null | awk '{print int($0) - 1}')

###############################################################################
# Format FASTQ files
###############################################################################

mkdir -p data/fastq

find . -type f |
    grep -e fastq -e fq |
    sed "s|_S.*||" |
    sort -u |
    while read -r line; do
        output="$(basename $line)"
        cat "$line"*R1* >data/fastq/"$output"_R1.fq.gz
        cat "$line"*R2* >data/fastq/"$output"_R2.fq.gz
    done

## FASTQC

mkdir -p reports/fastqc

for file in data/fastq/*; do
    fastqc -t "$threads" "$file" \
        --nogroup -o reports/fastqc
done

## MD5

md5sum data/fastq/*.fq.gz >reports/md5_fastq.txt
