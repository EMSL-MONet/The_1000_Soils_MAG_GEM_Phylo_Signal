# ============================================================
# 02_mantel_phylo_signal.R
#
# Tests whether predicted exchange of each compound is phylogenetically
# conserved, separately for Bacteria and Archaea.
#
# Per compound:
#   - D_phylo:    patristic distance between every pair of MAGs, from the
#                 IQ-TREE ultrafast-bootstrap consensus tree with GTDB
#                 reference tips pruned (computed once per domain).
#   - D_compound: |difference| in symmetric-min-max-normalized flux
#                 between the same MAG pairs.
#   - vegan::mantel(D_phylo, D_compound), Spearman, 999 permutations.
#     Positive r = closer relatives exchange the compound more similarly.
#
# Zero-variance compounds (same flux in every MAG) give r = NA and are
# reported as untested. P-values are Benjamini-Hochberg-adjusted across
# the tested compounds within each domain (p.adjust ignores NAs).
#
# Note: the eight trace ions (Mg, Co2+, K+, Mn2+, Cl-, Ca2+, Cu2+, Zn2+)
# have identical flux profiles, and Enterobactin/Fe-enterochelin have
# identical distance matrices; all are kept as separate tests, which
# makes the BH correction conservative.
#
# The seed is set once before looping over compounds, so permutation
# p-values depend on compound order in the input table.
#
# Runtime: ~25 min for Bacteria (497 MAGs), ~1 min for Archaea.
# ============================================================

source(file.path("R", "config.R"))
source(file.path("R", "utils.R"))
suppressPackageStartupMessages(library(vegan))

mantel_phylo_signal <- function(tree_file, fba_file, out_file, domain) {
  fba  <- readr::read_tsv(fba_file, show_col_types = FALSE)
  flux <- normalized_flux_matrix(fba)
  D_phylo <- mag_patristic_distances(tree_file, colnames(flux))

  message(sprintf("%s: %d MAGs, %d compounds", domain, ncol(flux), nrow(flux)))

  set.seed(PARAMS$seed)
  res <- dplyr::bind_rows(lapply(seq_len(nrow(flux)), function(i) {
    D_compound <- as.matrix(stats::dist(flux[i, ]))
    m <- vegan::mantel(D_phylo, D_compound,
                       method = PARAMS$mantel_method,
                       permutations = PARAMS$mantel_permutations)
    tibble::tibble(compound_id = fba$compound_id[i], name = fba$name[i],
                   mantel_r = unname(m$statistic), p_value = m$signif)
  }))

  res <- res %>%
    mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
    arrange(desc(mantel_r))

  message(sprintf("%s: %d tested, %d significant (P_adj < %.2f)",
                  domain, sum(!is.na(res$mantel_r)),
                  sum(res$p_adj < PARAMS$padj_cutoff, na.rm = TRUE), PARAMS$padj_cutoff))
  readr::write_tsv(res, out_file)
  invisible(res)
}

bacteria <- suppressWarnings(mantel_phylo_signal(
  INPUT$tree_bacteria, OUTPUT$fba_filtered_bacteria, OUTPUT$mantel_bacteria, "Bacteria"))
archaea <- suppressWarnings(mantel_phylo_signal(
  INPUT$tree_archaea, OUTPUT$fba_filtered_archaea, OUTPUT$mantel_archaea, "Archaea"))

# ---- Supplementary Table 2: tested bacterial compounds, ordered by r ----
bacteria %>%
  filter(!is.na(mantel_r)) %>%
  transmute(compound = clean_compound_name(name),
            r    = sprintf("%.3f", mantel_r),
            Padj = sprintf("%.3f", p_adj)) %>%
  readr::write_tsv(OUTPUT$supp_table_2)
