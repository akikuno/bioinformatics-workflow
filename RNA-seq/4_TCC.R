###############################################################################
# Input arguments
###############################################################################
file <- "count/count_gene_name.txt"
sample_name <- c("hetero-1", "hetero-2", "hetero-3", "ko-1", "ko-2", "ko-3")
group <- c(1, 1, 1, 2, 2, 2)
qval <- 0.05

###############################################################################
# Initialization
###############################################################################

options(repos = "http://cran.us.r-project.org")
if (!requireNamespace("pacman", quietly = T)) install.packages("pacman")
pacman::p_load(tidyverse, BiocManager, TCC, writexl)
system("mkdir -p degs")

###############################################################################
# Import data
###############################################################################

data <- read_tsv(file, skip = 1) %>% mutate(gene_id = toupper(Geneid))
colsubs <- seq(sample_name) + 6
colnames(data)[colsubs] <- sample_name

data <- data[rowSums(data[colsubs]) > 0, ]

###############################################################################
# TCC
###############################################################################

data_tcc <- data[, colsubs] %>% as.matrix()
rownames(data_tcc) <- toupper(data$Geneid)

tcc <- new("TCC", data_tcc, group) %>%
  filterLowCountGenes() %>%
  calcNormFactors(norm.method = "tmm", test.method = "edger", iteration = 3) %>%
  estimateDE(test.method = "edger", FDR = qval)

tmp_df <- data %>%
  select(gene_id, sample_name) %>%
  pivot_longer(-gene_id, names_to = "sample", values_to = "value") %>%
  mutate(group = str_remove(sample, "-[123]")) %>%
  group_by(gene_id, group) %>%
  mutate(mean = mean(value)) %>%
  select(gene_id, group, mean) %>%
  ungroup(group) %>%
  slice_max(mean, n = 1) %>%
  mutate(high_exp = group) %>%
  ungroup() %>%
  select(gene_id, high_exp) %>%
  distinct()

tcc_result <- getResult(tcc, sort = TRUE) %>%
  as_tibble() %>%
  mutate(estimatedDEG = if_else(q.value < qval, 1, 0)) %>%
  inner_join(tmp_df) %>%
  group_by(gene_id) %>%
  arrange(desc(estimatedDEG), q.value) %>%
  select(-rank)

sum(tcc_result$estimatedDEG == 1)

tcc_deg_hetero <- tcc_result %>% filter(estimatedDEG == 1, high_exp == "hetero")
tcc_deg_ko <- tcc_result %>% filter(estimatedDEG == 1, high_exp == "ko")

write_xlsx(list(all = tcc_result, hetero_up = tcc_deg_hetero, ko_up = tcc_deg_ko), "degs/degs.xlsx")
write_csv(tcc_result, "degs/degs.csv")