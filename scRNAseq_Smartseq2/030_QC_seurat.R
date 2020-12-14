#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#! Setting Options and Install Packages
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
options(repos = "http://cran.us.r-project.org")
if (!requireNamespace("pacman", quietly = T)) install.packages("pacman")
pacman::p_load(tidyverse, Seurat, patchwork, tidyseurat)

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#! load Data
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#==========================================================
#? Input
#==========================================================

data <- read_tsv("counts/counts.txt")

data_mat <- data[, -c(1:6)] %>% as.matrix
rownames(data_mat) <- data[, 1] %>% pull(Geneid)

data_meta <- data.frame(
    id = colnames(data_mat),
    celltype = colnames(data_mat) %>% str_remove("-.*"))

rownames(data_meta) <- data_meta$id
data_meta <- data_meta %>% select(-id)


data <- CreateSeuratObject(counts = data_mat, meta.data = data_meta)

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#! QC
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

data_qc <- data
data_qc <- PercentageFeatureSet(data_qc, "^MT-", col.name = "percent_mito")
data_qc <- PercentageFeatureSet(data_qc, "^RP[SL]", col.name = "percent_ribo")

feats <- c("nFeature_RNA","nCount_RNA","percent_mito","percent_ribo")
p_qc <- VlnPlot(data_qc, group.by= "celltype", features = feats, pt.size = 0.1,ncol = 4) + NoLegend()

(p_qc)

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#! UMAP
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
data_umap <- data

# Run the standard workflow for visualization and clustering
data_umap <- FindVariableFeatures(object = data_umap,
        selection.method = "vst", nfeatures = 2000, verbose = FALSE) %>%
    ScaleData(verbose = FALSE) %>%
    RunPCA(npcs = 30, verbose = FALSE) %>%
    RunUMAP(reduction = "pca", dims = 1:30)

p_umap <- DimPlot(object = data_umap, reduction = "umap", group.by = "celltype",
    label = TRUE, repel = TRUE) + NoLegend()

(p_umap)


#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#! PCA
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

data_pca <- data

# Run the standard workflow for visualization and clustering
data_pca <- FindVariableFeatures(object = data_pca,
        selection.method = "vst", nfeatures = 2000, verbose = FALSE) %>%
    ScaleData(verbose = FALSE) %>%
    RunPCA(npcs = 30, verbose = FALSE)

p_pca <- DimPlot(object = data_pca, reduction = "pca", group.by = "celltype",
    label = TRUE, repel = TRUE) + NoLegend()

(p_pca)


#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#! Outputs
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

ggsave((p_qc) / (p_umap + p_pca), filename = "results/qc_umap_pca.png", width = 10, height = 7)

data_qc <- subset(data_qc,
    subset = nFeature_RNA > 200 &
            percent_mito < 20)

data_qc <- NormalizeData(data_qc, normalization.method = "LogNormalize", scale.factor = 10000)

data_qc[["RNA"]]@data %>%
    as.data.frame() %>%
    rownames_to_column("gene_symbol") %>%
    as_tibble %>%
    write_csv("results/counts_normalized.csv")