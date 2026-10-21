# ============================================================
# utils.R
#
# Helpers shared across scripts.
# ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ape)
  library(treeio)
})

# Symmetric min-max scaling: x / max(|x|). Range [-1, 1]; preserves the
# sign (uptake vs. secretion) and true zeros. All-zero input returns zeros.
symmetric_minmax <- function(x) {
  m <- max(abs(x), na.rm = TRUE)
  if (m == 0) return(rep(0, length(x)))
  x / m
}

# Genome columns in an FBA table are headed "<phylum>||<MAG name>";
# any other column (compound_id, name, avg_uptake) is metadata.
genome_columns <- function(fba) {
  names(fba)[grepl("||", names(fba), fixed = TRUE)]
}

mag_names <- function(cols) sub("^.*\\|\\|", "", cols)

read_tree <- function(tree_file) {
  tree <- treeio::read.newick(tree_file)
  if (!inherits(tree, "phylo")) tree <- tree@phylo
  tree
}

# Patristic distance matrix among MAG tips only (reference tips pruned),
# rows/columns ordered to match `mags`.
mag_patristic_distances <- function(tree_file, mags) {
  tree <- read_tree(tree_file)
  missing <- setdiff(mags, tree$tip.label)
  if (length(missing) > 0) {
    stop("MAGs not found as tree tips: ", paste(missing, collapse = ", "))
  }
  d <- ape::cophenetic.phylo(ape::keep.tip(tree, mags))
  d[mags, mags]
}

# Flux matrix (compounds x MAGs) from a filtered FBA table, with each
# compound's row scaled by symmetric_minmax().
normalized_flux_matrix <- function(fba) {
  cols <- genome_columns(fba)
  m <- as.matrix(fba[, cols])
  out <- t(apply(m, 1, symmetric_minmax))
  dimnames(out) <- list(fba$compound_id, mag_names(cols))
  out
}

# Compound display name without the KBase compartment suffix ("[e0]", "[c0]").
clean_compound_name <- function(x) sub("\\s*\\[[a-z]0\\]$", "", x)
