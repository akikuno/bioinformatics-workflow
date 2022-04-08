#!/bin/sh

###############################################################################
# Define threads
###############################################################################

[ -z "$threads" ] && threads=$(getconf _NPROCESSORS_ONLN 2>/dev/null | awk '{print int($0) - 1}')
[ -z "$threads" ] && threads=1

###############################################################################
# Counting reads to genomic features by featureCounts
###############################################################################

mkdir -p reports

gtf="$(find data/mouse_genome/*gtf)"

featureCounts -T "$threads" -t exon -g gene_name -a "$gtf" \
  -o reports/raw_gene_counts_matrix.txt data/bam/*.bam

gzip reports/raw_gene_counts_matrix.txt

## MD5

md5sum reports/raw_gene_counts_matrix.txt.gz >reports/md5_counts.txt
