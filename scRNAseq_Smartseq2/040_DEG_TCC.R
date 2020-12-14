#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#! Setting Options and Install Packages
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
options(repos = "http://cran.us.r-project.org")
if (!require("pacman", quietly = TRUE)) install.packages("pacman")
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
pacman::p_load(tidyverse, TCC, Seurat, patchwork)

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#! load Data
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#==========================================================
#? Input
#==========================================================

data <- read_csv("counts/counts_normalized.csv")

data_mat <- as.matrix(data[, -1])
rownames(data_mat) <- pull(data[, 1])

group <- colnames(data)[-1] %>%
    str_remove_all("-.*") %>%
    as.factor() %>%
    as.integer()

tcc <- new("TCC", data_mat, group)
tcc <- calcNormFactors(tcc, norm.method = "tmm", test.method = "edger", iteration = 1)
tcc <- estimateDE(tcc, test.method = "edger", FDR = 0.1)
result <- getResult(tcc, sort = TRUE)

head(result)
sum(result$estimatedDEG == 1)