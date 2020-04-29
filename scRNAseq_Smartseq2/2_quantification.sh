#!/bin/sh

# ##################################################
# Define threads
# ##################################################
# Linux and similar...
[ -z "$threads" ] && threads=$(getconf _NPROCESSORS_ONLN 2>/dev/null | awk '{print int($0/2)}')
# FreeBSD and similar...
[ -z "$threads" ] && threads=$(getconf NPROCESSORS_ONLN | awk '{print int($0/2)}')
# Solaris and similar...
[ -z "$threads" ] && threads=$(ksh93 -c 'getconf NPROCESSORS_ONLN' | awk '{print int($0/2)}')
# Give up...
[ -z "$threads" ] && threads=1

# ##################################################
# featureCounts
# ##################################################

mkdir -p counts

bam_file=$(find . -name "*sorted.bam" -type f)

time featureCounts -T "$threads" -p -t exon -g gene_name \
    -a mouse_genome/Mus_musculus.GRCm38.99.gtf \
    -o counts/counts.txt \
    $(echo "$bam_file")

cat counts/counts.txt |
sed 1d |
sed -e "s#_Aligned.out_sorted.bam##g" \
    -e "s#./bam/##g" |
cat - > counts/counts_trimmed.txt