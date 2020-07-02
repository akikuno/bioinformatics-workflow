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
