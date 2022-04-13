###############################################################################
# Initialization
###############################################################################

options(repos = "http://cran.us.r-project.org")
if (!require("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(tidyverse, writexl, BiocManager, TCC)

###############################################################################
# Input arguments
###############################################################################

file <- "reports/raw_gene_counts_matrix.txt.gz"
samples <- c("cko", "ff")
target <- "cko"
nums <- c(4, 4)
qval <- 0.05


sample_name <- vector()
group <- vector()
for (i in seq_along(nums)) {
    sample_name <- append(sample_name, map2_chr(samples[i], seq(nums[i]), ~ str_c(.x, .y, sep = "-")))
    group <- append(group, rep(i, nums[i]))
}


###############################################################################
# Import data
###############################################################################

data <- read_tsv(file, skip = 1) %>% mutate(gene_id = toupper(Geneid))
colsubs <- 7:(sum(nums) + 6)
colnames(data)[colsubs] <- sample_name

data <- data[rowSums(data[colsubs]) > 0, ]

###############################################################################
# TCC
###############################################################################

data_tcc <- data[, colsubs] %>% as.matrix()
rownames(data_tcc) <- toupper(data$Geneid)

tcc <- new("TCC", data_tcc, group) %>%
    calcNormFactors(norm.method = "tmm", test.method = "edger", iteration = 3) %>%
    estimateDE(test.method = "edger", FDR = qval)

compare_expression <- data %>%
    select(gene_id, sample_name) %>%
    pivot_longer(-gene_id, names_to = "sample", values_to = "value") %>%
    mutate(genotype = str_remove(sample, "-[0-9]+")) %>%
    group_by(gene_id, genotype) %>%
    mutate(mean = mean(value)) %>%
    select(gene_id, genotype, mean) %>%
    ungroup(genotype) %>%
    distinct() %>%
    slice_max(mean, n = 1) %>%
    mutate(expression = case_when(
        genotype == target ~ str_c(target, "-high"),
        genotype != target ~ str_c(target, "-low")
    )) %>%
    ungroup() %>%
    select(gene_id, expression) %>%
    distinct()

tcc_result <- getResult(tcc, sort = TRUE) %>%
    as_tibble() %>%
    mutate(estimatedDEG = if_else(q.value < qval, 1, 0)) %>%
    inner_join(compare_expression) %>%
    group_by(gene_id) %>%
    arrange(desc(estimatedDEG), q.value) %>%
    select(-rank) %>%
    ungroup()

sum(tcc_result$estimatedDEG == 1)

###############################################################################
# Export results
###############################################################################

tcc_list <- tcc_result %>%
    filter(estimatedDEG == 1) %>%
    group_split(expression) %>%
    as.list()

for (i in seq_along(tcc_list)) {
    names(tcc_list)[i] <- tcc_list[[i]]$expression[1]
}
tcc_list <- append(list(all = tcc_result), tcc_list)

write_xlsx(
    tcc_list,
    "reports/degs.xlsx"
)

write_csv(tcc_result, "reports/degs.csv")