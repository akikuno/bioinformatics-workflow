###############################################################################
# Initialization
###############################################################################

options(repos = "http://cran.us.r-project.org")
if (!requireNamespace("pacman", quietly = T)) install.packages("pacman")
pacman::p_load(tidyverse, ggrepel)

###############################################################################
# Input arguments
###############################################################################

file <- "reports/degs.csv"

###############################################################################
# Import data
###############################################################################

tcc_result <- read_csv(file)

annotation <- c("ko-low", "ko-high")
legends <- c("KO DOWN", "KO UP")

###############################################################################
# Volucano and MA plots
###############################################################################

qval_min10 <- sort(tcc_result$q.value)[10]
mval_min10 <- sort(tcc_result$m.value)[10]
mval_max10 <- sort(tcc_result$m.value, decreasing = TRUE)[10]

tcc_volcano <-
    tcc_result %>%
    mutate(log10qval = -log10(q.value)) %>%
    mutate(label = case_when(
        estimatedDEG == 1 & expression == annotation[1] ~ legends[1],
        estimatedDEG == 1 & expression == annotation[2] ~ legends[2],
        TRUE ~ "NO"
    )) %>%
    mutate(annotate = case_when(
        estimatedDEG == 1 & str_detect(gene_id, TARGETS) ~ gene_id,
        TRUE ~ as.character(NA)
    ))


mycolors <- c("#5588BB", "#F06060", "#E5E5E5")
names(mycolors) <- c(legends[1], legends[2], "NO")


g_volcano <-
    ggplot(tcc_result, aes(x = m.value, y = log10qval, color = label, label = annotate)) +
    geom_point() +
    geom_text_repel() +
    scale_colour_manual(values = mycolors) +
    labs(x = "log2 Fold Change", y = "-log10 Q-value", colour = "") +
    theme_bw() +
    theme(
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 15),
        legend.text = element_text(size = 15),
        plot.background = element_rect(fill = "white")
    )


###############################################################################
# Export results
###############################################################################

ggsave("reports/volcano.png", g_volcano, width = 10, height = 5)
ggsave("reports/volcano.pdf", g_volcano, width = 10, height = 5)