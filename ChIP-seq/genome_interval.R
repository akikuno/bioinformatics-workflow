if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(tidyverse, janitor, BiocManager, IRange, fuzzyjoin)

klf4 <- tibble(
  id1 = 1:4,
  chromosome = c("chr1", "chr1", "chr2", "chr2"),
  start = c(100, 200, 300, 400),
  end = c(150, 250, 350, 450)
)

klf4 <- klf4 %>%
  mutate(center = start + (end - start) / 2)

oct4 <- tibble(
  id2 = 1:4,
  chromosome = c("chr1", "chr2", "chr2", "chr1"),
  start = c(140, 210, 400, 300),
  end = c(160, 240, 415, 320)
)

tmp_join <- genome_inner_join(klf4, oct4, by = c("chromosome", "start", "end"))

tmp_join %>%
  select(center, ends_with("y")) %>%
  janitor::clean_names()
rena()

# other functions:
genome_full_join(x1, x2, by = c("chromosome", "start", "end"))
genome_left_join(x1, x2, by = c("chromosome", "start", "end"))
genome_right_join(x1, x2, by = c("chromosome", "start", "end"))
genome_semi_join(x1, x2, by = c("chromosome", "start", "end"))
genome_anti_join(x1, x2, by = c("chromosome", "start", "end"))

# > # A tibble: 2 x 4
# >     id1 chromosome start   end
# >   <int> <chr>      <dbl> <dbl>
# > 1     2 chr1         200   250
# > 2     3 chr2         300   350