#!/usr/bin/env Rscript
## Arguments
input <- commandArgs(trailingOnly = TRUE)[1]
output <- commandArgs(trailingOnly = TRUE)[2]

## Options and Packages
options(repos = "https://cran.ism.ac.jp/")
options(scipen = 100)
if (!requireNamespace("pacman", quietly = T)) install.packages("pacman")
pacman::p_load(tidyverse)

# ++++++++++++++++++++++++++++++++++++++++++++++
## Read featureCounts data ----
# ++++++++++++++++++++++++++++++++++++++++++++++

counts <- read_tsv(input, skip = 1, col_types = cols())
col_name <- colnames(counts[, -c(1:5)]) %>%
    str_remove(".*/") %>%
    str_remove("_Aligned.*$") %>%
    str_c("TPM_", .)
# ++++++++++++++++++++++++++++++++++++++++++++++
## Calculate TPM ----
# ++++++++++++++++++++++++++++++++++++++++++++++

mat <- counts[, -c(1:5)] %>% as.matrix()
colnames(mat) <- col_name
rownames(mat) <- pull(counts[, 1])
rpk <- mat[, -1] / mat[, 1] * 1000
tpm <- apply(rpk, 2, function(x) x / sum(x) * 1000000)
colSums(tpm)

tpm <- data.frame(counts[, 1], tpm)
write_csv(x = data.frame(tpm), path = output)