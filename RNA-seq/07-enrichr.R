###############################################################################
# Initialization
###############################################################################

options(repos = "http://cran.us.r-project.org")
if (!requireNamespace("pacman", quietly = T)) install.packages("pacman")
pacman::p_load(tidyverse, enrichR)

###############################################################################
# Input arguments
###############################################################################

file <- "reports/degs.csv"

pathway <- c("WikiPathway_2021")

###############################################################################
# Import data
###############################################################################

tcc_result <- read_csv(file)

###############################################################################
# Enrichr
###############################################################################

dbs <- listEnrichrDbs() %>%
    as_tibble() %>%
    filter(str_detect(libraryName, pathway))

mcto_up <- tcc_result %>%
    filter(estimatedDEG == 1 & high_exp == "mcto") %>%
    pull(gene_id)
mcto_down <- tcc_result %>%
    filter(estimatedDEG == 1 & high_exp == "wt") %>%
    pull(gene_id)

enrichr_mcto_up <- enrichr(mcto_up, dbs)
enrichr_mcto_down <- enrichr(mcto_down, dbs)

df_mcto_up <- enrichr_mcto_up[[dbs$libraryName]] %>%
    as_tibble() %>%
    filter(Adjusted.P.value < 0.05) %>%
    head(10) %>%
    select(Term, Adjusted.P.value) %>%
    mutate(Term = str_remove(Term, " WP[0-9].*$")) %>%
    mutate(log10AdjPval = -log10(Adjusted.P.value))

df_mcto_down <- enrichr_mcto_down[[dbs$libraryName]] %>%
    as_tibble() %>%
    filter(Adjusted.P.value < 0.05) %>%
    head(10) %>%
    select(Term, Adjusted.P.value) %>%
    mutate(Term = str_remove(Term, " WP[0-9].*$")) %>%
    mutate(log10AdjPval = -log10(Adjusted.P.value))

###############################################################################
# Bar plot
###############################################################################

g_mcto_up <-
    ggplot(df_mcto_up, aes(x = fct_reorder(Term, log10AdjPval), y = log10AdjPval)) +
    geom_col() +
    scale_y_reverse() +
    coord_flip() +
    labs(x = "-log10 Adjusted P-value", y = "", colour = "") +
    theme_bw() +
    theme(
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 15),
        legend.text = element_text(size = 15)
    )

g_mcto_down <-
    ggplot(df_mcto_down, aes(x = fct_reorder(Term, log10AdjPval), y = log10AdjPval)) +
    geom_col() +
    scale_y_reverse() +
    coord_flip() +
    labs(x = "-log10 Adjusted P-value", y = "", colour = "") +
    theme_bw() +
    theme(
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 15),
        legend.text = element_text(size = 15)
    )
###############################################################################
# Export results
###############################################################################

ggsave("reports/enrchr_mcto_up.png", g_mcto_up, width = 13, height = 8)
ggsave("reports/enrich_mcto_up.pdf", g_mcto_up, width = 13, height = 8)
ggsave("reports/enrchr_mcto_down.png", g_mcto_down, width = 13, height = 8)
ggsave("reports/enrich_mcto_down.pdf", g_mcto_down, width = 13, height = 8)