# ============================================================
# 05_figure2_mantel_scatter.R
#
# Figure 2 and Supplementary Figure 1: the two matrices compared by each
# whole-tree Mantel test, flattened to one value pair per MAG pair
# (x = patristic distance, y = |difference| in normalized flux), shown
# as 2D bin densities (40 x 40 bins, log10 count) with an OLS regression
# line and 95% confidence band for visual reference. Each panel is
# annotated with the Spearman Mantel r and BH-adjusted P from script 02.
# Panels are ordered by decreasing r.
#
#   Figure 2:              compounds with P_adj < 0.05
#   Supplementary Figure 1: all other tested compounds
# ============================================================

source(file.path("R", "config.R"))
source(file.path("R", "utils.R"))
suppressPackageStartupMessages(library(ggplot2))

fba    <- readr::read_tsv(OUTPUT$fba_filtered_bacteria, show_col_types = FALSE)
flux   <- normalized_flux_matrix(fba)
mantel <- readr::read_tsv(OUTPUT$mantel_bacteria, show_col_types = FALSE) %>%
  filter(!is.na(mantel_r)) %>%
  arrange(desc(mantel_r))

D_phylo <- mag_patristic_distances(INPUT$tree_bacteria, colnames(flux))
upper   <- upper.tri(D_phylo)

scatter_data <- function(ids) {
  dplyr::bind_rows(lapply(ids, function(id) {
    D_compound <- as.matrix(stats::dist(flux[id, ]))
    data.frame(compound_id = id, D_phylo = D_phylo[upper], D_compound = D_compound[upper])
  })) %>%
    left_join(mantel, by = "compound_id") %>%
    mutate(panel = sprintf("%s\nr = %.3f (Padj = %.3f)",
                           clean_compound_name(name), mantel_r, p_adj),
           panel = factor(panel, levels = unique(panel[order(-mantel_r)])))
}

scatter_plot <- function(ids, ncol) {
  ggplot(scatter_data(ids), aes(x = D_phylo, y = D_compound)) +
    geom_bin2d(bins = 40) +
    geom_smooth(method = "lm", formula = y ~ x, color = "red", se = TRUE, linewidth = 0.7) +
    scale_fill_viridis_c(trans = "log10", name = "MAG pairs\n(count)") +
    facet_wrap(~panel, ncol = ncol, scales = "free_y") +
    labs(x = "Patristic (phylogenetic) distance between MAG pair",
         y = "Normalized flux difference between MAG pair") +
    theme_minimal(base_size = 11) +
    theme(strip.text = element_text(size = 8, face = "bold"),
          axis.text  = element_text(size = 6))
}

sig   <- mantel$compound_id[mantel$p_adj <  PARAMS$padj_cutoff]
other <- mantel$compound_id[mantel$p_adj >= PARAMS$padj_cutoff]
message(sprintf("Figure 2: %d panels; Supplementary Figure 1: %d panels", length(sig), length(other)))

ggsave(file.path(FIG_DIR, "Figure_2_mantel_scatter_significant.pdf"),
       scatter_plot(sig, ncol = 4), width = 12, height = 3.2 * ceiling(length(sig) / 4))
ggsave(file.path(FIG_DIR, "Supplementary_Figure_1_mantel_scatter_other.pdf"),
       scatter_plot(other, ncol = 8), width = 22, height = 2.8 * ceiling(length(other) / 8))
