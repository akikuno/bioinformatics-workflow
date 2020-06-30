## Options and Packages
options(repos = "https://cran.ism.ac.jp/")
if (!requireNamespace("pacman", quietly = T)) install.packages("pacman")
if (!requireNamespace("BiocManager", quietly = T)) install.packages("BiocManager")
pacman::p_load(tidyverse, HiCBricks)

Bintable_path <- system.file(file.path("extdata",
"Bintable_100kb.bins"), package = "HiCBricks")

BrickContainer_dir <- file.path(tempdir(), "HiCBricks_vignette_test")
My_BrickContainer <- load_BrickContainer(project_dir = BrickContainer_dir)
Example_dataset_dir <- system.file("extdata", package = "HiCBricks")

Chromosomes <- c("chr2L", "chr3L", "chr3R", "chrX")
for (chr in Chromosomes) {
    Matrix_file <- file.path(Example_dataset_dir,
        paste(paste("Sexton2012_yaffetanay_CisTrans_100000_corrected", 
          chr, sep = "_"), "txt.gz", sep = "."))
    Brick_load_matrix(Brick = My_BrickContainer,
        chr1 = chr,
        chr2 = chr,
        resolution = 100000,
        matrix_file = Matrix_file,
        delim = " ",
        remove_prior = TRUE)
}