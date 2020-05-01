#!/usr/bin/env Rscript
## Arguments
input = commandArgs(trailingOnly=TRUE)[1]
output = commandArgs(trailingOnly=TRUE)[2]

## Options and Packages
options(repos = "https://cran.ism.ac.jp/")
options(scipen = 100)
if (!requireNamespace("pacman", quietly = T)) install.packages("pacman")
pacman::p_load(tidyverse)

# ++++++++++++++++++++++++++++++++++++++++++++++
## Read featureCounts data ----
# ++++++++++++++++++++++++++++++++++++++++++++++

counts <- read_tsv(input, skip = 1)

# ++++++++++++++++++++++++++++++++++++++++++++++
## Calculate TPM ----
# ++++++++++++++++++++++++++++++++++++++++++++++

df <- counts[,-c(1:5)]
colnames(df) <- str_c("TPM_", colnames(df))
rownames(df) <- counts[, 1]
rpk <- df[,-1]/df[,1] * 1000
tpm <- apply(rpk, 2, function(x) x/sum(x) * 1000000)
colSums(tpm)

tpm <- data.frame(counts[,1], tpm)
write_csv(x = data.frame(tpm), path = output)
