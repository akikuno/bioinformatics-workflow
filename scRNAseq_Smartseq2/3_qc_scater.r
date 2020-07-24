#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#! Setting Options and Install Packages
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
options(repos = "http://cran.us.r-project.org")
if (!requireNamespace("pacman", quietly = T)) install.packages("pacman")
pacman::p_load(tidyverse, scRNAseq, scater, Seurat, patchwork)


#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#! load Data
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#==========================================================
#? TEST Input
#==========================================================

sce.416b <- LunSpikeInData(which="416b")
sce.416b$block <- factor(sce.416b$block)

location <- rowRanges(sce.416b)
is.mito <- any(seqnames(location)=="MT")

library(scater)
df <- perCellQCMetrics(sce.416b, subsets=list(Mito=is.mito))
df


data_sce <- as.SingleCellExperiment(data)
altExpNames(data_sce) <- "ERCC"