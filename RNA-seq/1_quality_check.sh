#!/bin/sh

# ##################################################
# Define threads
# ##################################################

[ -z "$threads" ] && threads=$(getconf _NPROCESSORS_ONLN 2>/dev/null | awk '{print int($0/2)}')

# ##################################################
# Format FASTQ files
# ##################################################

mkdir -p fastq
find . -type f |
    grep -e fastq -e fq |
    sed "s|_S.*||" |
    sort -u |
    while read -r line; do
        output="$(basename $line)"
        cat "$line"*R1* >fastq/"$output"_R1.fq.gz
        cat "$line"*R2* >fastq/"$output"_R2.fq.gz
    done

## FASTQC
mkdir -p ./fastqc
for file in ./fastq/*; do
    fastqc -t "$threads" "$file" \
        --nogroup -o ./fastqc
done

# # ##################################################
# # fastp trimming
# # ##################################################

# R1=$(ls ./fastq/*R1*)
# R2=$(ls ./fastq/*R2*)

# num=$(find ./fastq/*R1* -type f | awk '{print NR}')
# #awk '{if (NR%2!=0) a++; print a }'

# mkdir -p fastq_trim

# for i in $(echo $num); do
#     fw=$(echo $R1 | cut -d " " -f $i)
#     rv=$(echo $R2 | cut -d " " -f $i)
#     echo $fw
#     #
#     out_fw=$(echo ${fw} | sed -e "s#.*/#fastq_trim/#g" -e "s/R1/R1_trim/g")
#     out_rv=$(echo ${rv} | sed -e "s#.*/#fastq_trim/#g" -e "s/R2/R2_trim/g")
#     report=$(echo ${out_fw} | sed "s/_R1*//g")
#     #
#     time fastp -i ${fw} -I ${rv} \
#         -o ${out_fw} -O ${out_rv} \
#         -f 5 -t 1 \
#         -h ${report}.html -j ${report}.json \
#         --thread ${threads}
# done

# ## FASTQC
# mkdir -p ./fastqc_trim_report
# for file in ./fastq_trim/*.gz; do
#     fastqc -t "$threads" "$file" \
#         --nogroup -o ./fastqc_trim_report
# done
