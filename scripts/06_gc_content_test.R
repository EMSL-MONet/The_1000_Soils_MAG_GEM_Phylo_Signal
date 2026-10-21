# ============================================================
# 06_gc_content_test.R
#
# Compares genomic GC content (CheckM2) between bacterial MAGs that
# combine branched-chain amino acid (BCAA) import with secretion of
# fermentation products and all other bacterial MAGs.
#
# Pattern definition (from the filtered FBA table):
#   uptake (> 0) of L-leucine and/or L-isoleucine, AND
#   secretion (< 0) of both pyruvate and ethanol.
#
# Two-sided Wilcoxon rank-sum (Mann-Whitney U) test.
# Expected: 36 pattern MAGs vs. 461 others, W = 7484, P = 0.33.
# ============================================================

source(file.path("R", "config.R"))
source(file.path("R", "utils.R"))

fba  <- readr::read_tsv(OUTPUT$fba_filtered_bacteria, show_col_types = FALSE)
cols <- genome_columns(fba)

flux <- function(compound) {
  setNames(unlist(fba[fba$name == compound, cols]), mag_names(cols))
}

pattern <- (flux("L-Leucine [e0]") > 0 | flux("L-Isoleucine [e0]") > 0) &
           flux("Pyruvate [e0]") < 0 & flux("Ethanol [e0]") < 0

checkm2 <- readr::read_tsv(INPUT$checkm2, show_col_types = FALSE)
gc <- setNames(checkm2$GC_Content, checkm2$Name)[names(pattern)]

test <- wilcox.test(gc[pattern], gc[!pattern], alternative = "two.sided", exact = FALSE)

result <- tibble::tibble(
  n_pattern       = sum(pattern),
  n_other         = sum(!pattern),
  gc_median_pattern = median(gc[pattern]),
  gc_min_pattern    = min(gc[pattern]),
  gc_max_pattern    = max(gc[pattern]),
  gc_median_other   = median(gc[!pattern]),
  W               = unname(test$statistic),
  p_value         = test$p.value
)
print(as.data.frame(result))
readr::write_tsv(result, OUTPUT$gc_test)
