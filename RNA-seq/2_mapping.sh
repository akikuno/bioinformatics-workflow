#!/bin/sh

###############################################################################
# Define threads
###############################################################################

[ -z "$threads" ] && threads=$(getconf _NPROCESSORS_ONLN 2>/dev/null | awk '{print int($0) - 1}')
[ -z "$threads" ] && threads=1

###############################################################################
# Genome indexing by Subread
###############################################################################

mkdir -p data/mouse_index/subread

time subread-buildindex \
  -o data/mouse_index/subread/subread \
  data/mouse_genome/Mus_musculus.GRCm38.dna.primary_assembly.fa

###############################################################################
# Mapping to genome by Subread
###############################################################################

mkdir -p data/bam

find data/fastq/ -type f |
  sort |
  awk 'NR%2==1 {printf $0" "; next}1' |
  while read -r R1 R2; do
    out_prefix="$(basename ${R1%%_R*})"
    #
    subread-align -t 0 \
      -T "${threads:-1}" \
      -i data/mouse_index/subread/subread \
      -r "$R1" -R "$R2" \
      -o data/bam/tmp.bam
    #
    samtools sort -@ "$threads" data/bam/tmp.bam >data/bam/"$out_prefix".bam
    samtools index -@ "$threads" data/bam/"$out_prefix".bam
    samtools stats -@ "$threads" data/bam/"$out_prefix".bam >data/bam/"$out_prefix"_stats
    rm data/bam/tmp.bam*
  done

###############################################################################
# BigWig files to visualize by IGV
###############################################################################

mkdir -p data/bigwig

for bam in bam/*bam; do
  out_f=$(echo "$bam" |
    sed "s|bam/|bigwig/|" |
    sed "s/.bam$/_bin5_cpm.bw/g")

  bamCoverage -b "$bam" -o "$out_f" -p "$threads" \
    --binSize 5 --normalizeUsing CPM
done

###############################################################################
# Multiqc
###############################################################################

multiqc .

# # ----------------------------------------------
# # past code
# # ----------------------------------------------

# # ##################################################
# # STAR index
# # ##################################################

# mkdir -p data/mouse_index/STAR

# STAR \
# --runMode genomeGenerate \
# --genomeDir data/mouse_index/STAR \
# --genomeFastaFiles data/mouse_genome/Mus_musculus.GRCm38.dna.primary_assembly.fa \
# --sjdbGTFfile data/mouse_genome/Mus_musculus.GRCm38.99.gtf \
# --runThreadN "$threads"

# # ##################################################
# # STAR mapping
# # ##################################################

# R1=$(ls ./fastq/*R1*.gz)
# R2=$(ls ./fastq/*R2*.gz)
# num=$(find ./fastq/*R1*.gz -type f | awk '{print NR}')

# mkdir -p bam

# for i in $(echo $num) ; do
#     fw=$(echo $R1 | cut -d " " -f $i)
#     rv=$(echo $R2 | cut -d " " -f $i)
#     out_f=$(echo "$fw" |
#     sed -e "s#.*/#bam/#g" -e "s/_R1.*//g")
#     #
#     time STAR --runThreadN "$threads" \
#     --genomeDir data/mouse_index/STAR \
#     --readFilesIn "$fw" "$rv" \
#     --readFilesCommand gunzip -c \
#     --outSAMtype BAM Unsorted \
#     --outFileNamePrefix "$out_f"_
# done

# # ##################################################
# # Sort and index
# # ##################################################

# for bam in ./bam/*bam ; do
#     out_f=$(echo $bam | sed "s/\.bam/_sorted.bam/g")
#     samtools sort -@ "$threads" $bam -o "$out_f"
#     samtools index -@ "$threads" "$out_f"
#     rm $bam
# done

# # ##################################################
# # Quantification by RSEM
# # ##################################################

# # Index
# mkdir -p data/mouse_index/RSEM/index
# rsem-prepare-reference \
# --star \
# --gtf data/mouse_genome/Mus_musculus.GRCm38.99.gtf \
# --num-threads "$threads" \
# data/mouse_genome/Mus_musculus.GRCm38.dna.primary_assembly.fa \
# data/mouse_index/RSEM/index

# # Quantification
# fastqs="fastq fastq_trim"
# bams="bam bam_trim"
# mkdir -p ${bams}

# for ijk in 1 2; do
#     bam=$(echo ${bams} | cut -d " " -f ${ijk})
#     fastq=$(echo ${fastqs} | cut -d " " -f ${ijk})
#     #
#     R1=$(ls ./${fastq}/*R1*.gz)
#     R2=$(ls ./${fastq}/*R2*.gz)
#     num=$(find ./${fastq}/*R1*.gz -type f | awk '{print NR}')
#     #
#     for i in $(echo ${num}) ; do
#         fw=$(echo ${R1} | cut -d " " -f ${i})
#         rv=$(echo ${R2} | cut -d " " -f ${i})
#         out_f=$(echo "${fw}" |
#         sed -e "s/${fastq}/${bam}/g" -e "s/_R1*.gz//g")
#         echo "${out_f} is now processing..."
#         #
#         { gzip -dc ${fw} > /tmp/R1* & } 2>/dev/null
#         { gzip -dc ${rv} > /tmp/R2* & } 2>/dev/null
#         wait 2>/dev/null
#         time rsem-calculate-expression \
#         -p ${threads} \
#         --paired-end \
#         --star \
#         --output-genome-bam \
#         /tmp/R1* /tmp/R2* \
#         data/mouse_index/RSEM/index \
#         ${out_f}
#     done
# done
# multiqc .
