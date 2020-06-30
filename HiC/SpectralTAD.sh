#!/bin/bash


###########################################################
#! Conda environment
###########################################################

# conda create -y -n r-env r-base=3.6.3 r-essentials=3.6
# conda install -y -n r-env -c bioconda bioconductor-spectraltad ucsc-liftover
conda activate r-env

###########################################################
#! SpectralTAD
###########################################################

cat << EOF > SpectralTAD.R
# SpectralTAD: https://bit.ly/3dOJtZE
## Options and Packages
options(repos = "https://cran.ism.ac.jp/")
if (!requireNamespace("pacman", quietly = T)) install.packages("pacman")
pacman::p_load(tidyverse, SpectralTAD)

args <- commandArgs(trailingOnly=TRUE)
input <- args[1]
output <- args[2]
df <- read_csv(input)
mat <- df
mat <- mat[,-1] %>% as.matrix
rownames(mat) <-  df[,1] %>% pull
mat[is.nan(mat)] <- 0.0001

result = SpectralTAD(mat, chr = "chrX", resolution = 40000, qual_filter = FALSE, z_clust = FALSE)
write_tsv(as.data.frame(result), output)
EOF
chmod +x SpectralTAD.R

inputs=$(
    find /mnt/d/nisimura-sensei_xi/Giorgetti_HiC/* |
    grep 40kb |
    grep iced-snpMasked |
    grep chrX |
    grep matrix.gz
)

# outputs=$(
#     find /mnt/d/nisimura-sensei_xi/Giorgetti_HiC/* |
#     grep 40kb |
#     grep iced-snpMasked |
#     grep chrX |
#     grep matrix.gz |
#     sed "s#.*/##g" |
#     sed "s/EHSNP-//g" |
#     sed "s/__mm9-cast-129s1__genome__C-40000-iced.masked__/-/g" |
#     sed "s/.matrix.gz/.bg/g"
# )


format()(
    cat - |
    awk 'NR==1{for(i=1;i<=NF; i++) sub(/.*-/, "", $i)}1' |
    sed "s/^.*-//g" |
    sed "s/[ |\t]/,/g"
)

echo "$inputs" |
while read -r input; do
    output=$(
        echo "$input" |
        sed "s#.*/##g" |
        sed "s/EHSNP-//g" |
        sed "s/__mm9-cast-129s1__genome__C-40000-iced.masked__/-/g" |
        sed "s/.matrix.gz/.bg/g"
    )
    echo "$output is generating..."
    gzip -dc "$input" |
    format > test.csv
    Rscript SpectralTAD.R test.csv "$output"
done
rm test.csv

# awk -v inputs="$inputs" -v outputs="$outputs" \
# 'BEGIN{split(inputs, inputs_, " ")
#     split(outputs, outputs_, " ")
#     for(i in inputs_) {
#         print "echo ", outputs_[i], "is now generating..."
#         print "Rscript SpectralTAD.R",inputs_[i], outputs_[i]
#         }
#     }' |
# sh -

###########################################################
#! LIFTOVER
###########################################################

wget http://hgdownload.cse.ucsc.edu/goldenpath/mm9/liftOver/mm9ToMm10.over.chain.gz

for bg in *bg; do
output=$(echo "$bg" | sed "s/.bg/_mm10.bg/g")
echo $output
cat "${bg}" | sed 1d > test.bg
liftOver test.bg mm9ToMm10.over.chain.gz "${output}" unlifted
done
rm test.bg unlifted