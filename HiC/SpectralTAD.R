# SpectralTAD: https://bit.ly/3dOJtZE

# Shellscript -----------------------------------
# conda create -y -n r-env r-base=3.6.3 r-essentials=3.6
# conda install -y -n r-env -c bioconda bioconductor-spectraltad
# conda activatre r-env
# R
# -----------------------------------------------

## Options and Packages
options(repos = "https://cran.ism.ac.jp/")
if (!requireNamespace("pacman", quietly = T)) install.packages("pacman")
pacman::p_load(tidyverse, SpectralTAD)

setwd("/mnt/d/nisimura-sensei_xi/Giorgetti_HiC/matrix__cis__txt__NPC/40kb/iced-snpMasked")
mat_129 <- list.files() %>% str_subset("chrX-129S1")
mat_cast <- list.files() %>% str_subset("chrX-cast")

mat_129 <- read_tsv(mat_129, col_names = F)
mat_cast <- read_tsv(mat_cast, col_names = F)

