#!/bin/bash

mkdir -p data/fastq reports

find . -type f |
    grep -e fastq -e fq |
    sed "s|_S.*||" |
    sort -u |
    while read -r line; do
        output="$(basename $line)"
        cat "$line"*R1* >data/fastq/"$output"_R1.fq.gz
        cat "$line"*R2* >data/fastq/"$output"_R2.fq.gz
    done

## MD5

md5sum data/fastq/*.fq.gz >reports/md5_fastq.txt
