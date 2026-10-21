# ============================================================
# run_all.R
#
# Runs the full analysis in order. From the repository root:
#   Rscript run_all.R
# Total runtime ~30 min (dominated by the bacterial Mantel tests).
# ============================================================

scripts <- c(
  "scripts/01_filter_fba_tables.R",
  "scripts/02_mantel_phylo_signal.R",
  "scripts/03_build_tip_metadata.R",
  "scripts/04_figure1_tree_with_rings.R",
  "scripts/05_figure2_mantel_scatter.R",
  "scripts/06_gc_content_test.R"
)

for (s in scripts) {
  message("\n==== ", s, " ====")
  source(s, local = new.env())
}

writeLines(capture.output(sessionInfo()), file.path("results", "sessionInfo.txt"))
