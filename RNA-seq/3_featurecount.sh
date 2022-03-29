#!/bin/sh

###############################################################################
# Counting reads to genomic features by featureCounts
###############################################################################

mkdir -p reports

gtf="$(find data/mouse_genome/*gtf)"

featureCounts -t exon -g gene_name -a "$gtf" \
  -o reports/count_gene_name.txt data/bam/*.bam

gzip -c reports/raw_gene_counts_matrix.txt >reports/raw_gene_counts_matrix.txt.gz

## MD5

md5sum reports/raw_gene_counts_matrix.txt.gz >reports/md5_counts.txt
