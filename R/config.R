# ============================================================
# config.R
#
# Input/output paths and analysis parameters used by analysis.Rmd.
# All paths are relative to the repository root; render
# analysis.Rmd from that folder (see README.md).
# ============================================================

if (!file.exists(file.path("R", "config.R"))) {
  stop("Render analysis.Rmd from the repository root (setwd() to it first).")
}

# ---- Inputs (data/) ----
INPUT <- list(
  # Compiled KBase FBA exchange fluxes: one row per compound, one column
  # per MAG, headed "<phylum>||<MAG name>". Positive = uptake,
  # negative = secretion, NA = compound absent from that MAG's model.
  fba_bacteria = file.path("data", "fba_uptake_bacteria.tsv"),
  fba_archaea  = file.path("data", "fba_uptake_archaea.tsv"),
  # IQ-TREE ultrafast-bootstrap consensus trees (MAGs + GTDB references)
  tree_bacteria = file.path("data", "bac120.contree"),
  tree_archaea  = file.path("data", "ar53.contree"),
  # CheckM2 quality report for the 551 dereplicated MAGs
  checkm2 = file.path("data", "checkm2_551_MAGs.tsv"),
  # Supplementary Table 2 as published, used to check the reproduction
  published_supp_table_2 = file.path("data", "published_supplementary_table_2.tsv")
)

# ---- Outputs (results/) ----
OUT_DIR <- "results"
FIG_DIR <- file.path(OUT_DIR, "figures")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

OUTPUT <- list(
  fba_filtered_bacteria = file.path(OUT_DIR, "fba_filtered_bacteria.tsv"),
  fba_filtered_archaea  = file.path(OUT_DIR, "fba_filtered_archaea.tsv"),
  mantel_bacteria       = file.path(OUT_DIR, "mantel_bacteria.tsv"),
  mantel_archaea        = file.path(OUT_DIR, "mantel_archaea.tsv"),
  supp_table_2          = file.path(OUT_DIR, "Supplementary_Table_2_bacteria_mantel_r_Padj.tsv"),
  tip_metadata_bacteria = file.path(OUT_DIR, "tip_metadata_bacteria.tsv"),
  gc_test               = file.path(OUT_DIR, "gc_content_exchange_pattern_test.tsv"),
  figure_1              = file.path(FIG_DIR, "Figure_1_bacteria_tree_with_rings.pdf"),
  figure_2              = file.path(FIG_DIR, "Figure_2_mantel_scatter_significant.pdf"),
  supp_figure_1         = file.path(FIG_DIR, "Supplementary_Figure_1_mantel_scatter_other.pdf")
)

# ---- Parameters ----
PARAMS <- list(
  missing_fraction_cutoff = 0.333, # drop compounds missing (NA) in >= this fraction of MAGs
  mantel_method           = "spearman",
  mantel_permutations     = 999,
  seed                    = 1,
  padj_cutoff             = 0.05,
  root_phylum_bacteria    = "p__Patescibacteriota"
)
