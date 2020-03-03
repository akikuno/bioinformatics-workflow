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
# Subread index
# ##################################################
# --------------------------------------------------
output_dir="mouse_index/subread"
# --------------------------------------------------
mkdir -p ${output_dir}
subread-buildindex \
    -o ${output_dir}/subread \
    mouse_genome/Mus_musculus.GRCm38.dna.primary_assembly.fa

# ##################################################
# Subread mapping
# ##################################################
# --------------------------------------------------
input_label="trim"
output_dir="bam_trim"
# --------------------------------------------------
mkdir -p ${output_dir}

R1=$(ls ./fastq/*R1*.gz)
R2=$(ls ./fastq/*R2*.gz)
num=$(find ./fastq/*R1*.gz -type f | awk '{print NR}')

for i in $(echo $num) ; do
    fw=$(echo $R1 | cut -d " " -f $i)
    rv=$(echo $R2 | cut -d " " -f $i)
    out_f=$(echo "$fw" |
    sed -e "s#.*/#bam/#g" -e "s/_R1.*/${input_label}/g")
    #
    echo "${out_f} is now processing..."
    subread-align -t 0 -T "$threads" -d 50 -D 600 -i mouse_index/subread/subread \
    -r ${fw} -R ${rv} \
    -o tmp_${out_f}.bam
    #
    samtools sort -@ "$threads" tmp.bam > "$out_f".bam
    samtools index -@ "$threads" "$out_f".bam
    samtools stats -@ "$threads" "$out_f".bam > ${out_f}_stats
    rm tmp_${out_f}.bam
done

# ##################################################
# featureCounts
# ##################################################
# --------------------------------------------------
bam_dir="bam_trim"
output="counts_gene_id.txt"
gtf="mouse_genome/Mus_musculus.GRCm38.99.gtf"
# --------------------------------------------------

featureCounts -t exon -g gene_id -a ${gtf} \
    -o ${output} ${bam_dir}/*.bam

# ##################################################
# BigWig files to visualize by IGV
# ##################################################
# --------------------------------------------------
bam_dir="bam_trim"
output_dir="bigwig"
# --------------------------------------------------

mkdir -p ${output_dir}

for bam in ./${bam_dir}/*bam ; do
    out_f=$(echo "$bam" |
        sed -e "s#${bam_dir}/#${bigwig}/#g" \
            -e "s/.bam$/_bin5_cpm.bw/g")
    bamCoverage -b $bam -o "$out_f" \
    -p "$threads" \
    --binSize 5 \
    --normalizeUsing CPM
done

# ##################################################
# Multiqc
# ##################################################
multiqc .


# # ----------------------------------------------
# # past code
# # ----------------------------------------------

# # ##################################################
# # STAR index
# # ##################################################

# mkdir -p mouse_index/STAR

# STAR \
# --runMode genomeGenerate \
# --genomeDir mouse_index/STAR \
# --genomeFastaFiles mouse_genome/Mus_musculus.GRCm38.dna.primary_assembly.fa \
# --sjdbGTFfile mouse_genome/Mus_musculus.GRCm38.99.gtf \
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
#     --genomeDir mouse_index/STAR \
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
# mkdir -p mouse_index/RSEM/index
# rsem-prepare-reference \
# --star \
# --gtf mouse_genome/Mus_musculus.GRCm38.99.gtf \
# --num-threads "$threads" \
# mouse_genome/Mus_musculus.GRCm38.dna.primary_assembly.fa \
# mouse_index/RSEM/index

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
#         mouse_index/RSEM/index \
#         ${out_f}
#     done
# done
# multiqc .
