#!/bin/sh

###############################################################################
# Counting reads to genomic features by featureCounts
###############################################################################

mkdir -p count

gtf="$(find mouse_genome/*gtf)"

featureCounts -t exon -g gene_name -a "$gtf" \
  -o count/count_gene_name.txt bam/*.bam

gzip -c count/count_gene_name.txt >count/count_gene_name.txt.gz
