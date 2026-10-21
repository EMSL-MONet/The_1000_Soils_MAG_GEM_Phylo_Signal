# ============================================================
# analysis.R
#
# Functions that reproduce the published analyses. Called from
# analysis.Rmd; each function returns its result and, where noted,
# also writes it to results/.
# ============================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggtree)
  library(ggtreeExtra)
  library(ggnewscale)
  library(tidyr)
  library(scales)
  library(vegan)
})

# ---- Step 1: filter the compiled FBA tables -------------------------------
#
#   1. Drop MAGs with no FBA output (every compound NA). For Bacteria this
#      removes PPRH-TOP_bin.5.o, whose draft model could not be gapfilled
#      to produce biomass.
#   2. Drop compounds missing (NA) in >= one-third of the remaining MAGs
#      (threshold = ceiling(0.333 x number of MAGs)).
#   3. Set the remaining NAs to zero (no exchange by that MAG's model).
#   4. Sort compounds by mean flux, highest first; ties keep compound-ID order.
filter_fba_table <- function(in_file, out_file = NULL) {
  fba  <- readr::read_tsv(in_file, show_col_types = FALSE)
  cols <- genome_columns(fba)

  no_model <- cols[colSums(!is.na(fba[, cols])) == 0]
  cols     <- setdiff(cols, no_model)

  threshold <- ceiling(length(cols) * PARAMS$missing_fraction_cutoff)
  n_missing <- rowSums(is.na(fba[, cols]))
  kept <- fba[n_missing < threshold, c("compound_id", "name", cols)]
  kept[cols] <- lapply(kept[cols], function(x) replace(x, is.na(x), 0))
  kept$avg_uptake <- rowMeans(kept[, cols])
  kept <- kept[order(-kept$avg_uptake), ]

  if (!is.null(out_file)) readr::write_tsv(kept, out_file)
  attr(kept, "summary") <- list(
    mags_dropped = mag_names(no_model), n_mags = length(cols),
    n_compounds_in = nrow(fba), n_compounds_kept = nrow(kept), threshold = threshold
  )
  kept
}

# ---- Step 2: Mantel test of phylogenetic signal ---------------------------
#
# Per compound: Spearman Mantel test (999 permutations) between pairwise
# patristic distance (consensus tree, reference tips pruned) and pairwise
# |difference| in symmetric-min-max-normalized flux. Zero-variance compounds
# give r = NA (untested). P-values are BH-adjusted across tested compounds.
#
# The seed is set once before looping over compounds, so p-values depend on
# the compound order of `fba` (kept fixed by filter_fba_table()).
mantel_phylo_signal <- function(tree_file, fba, out_file = NULL) {
  flux    <- normalized_flux_matrix(fba)
  D_phylo <- mag_patristic_distances(tree_file, colnames(flux))

  set.seed(PARAMS$seed)
  res <- dplyr::bind_rows(lapply(seq_len(nrow(flux)), function(i) {
    D_compound <- as.matrix(stats::dist(flux[i, ]))
    m <- suppressWarnings(vegan::mantel(D_phylo, D_compound,
                                        method = PARAMS$mantel_method,
                                        permutations = PARAMS$mantel_permutations))
    tibble::tibble(compound_id = fba$compound_id[i], name = fba$name[i],
                   mantel_r = unname(m$statistic), p_value = m$signif)
  }))

  res <- res %>%
    mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
    arrange(desc(mantel_r))

  if (!is.null(out_file)) readr::write_tsv(res, out_file)
  res
}

# Supplementary Table 2: tested bacterial compounds, ordered by r.
supplementary_table_2 <- function(mantel, out_file = NULL) {
  tab <- mantel %>%
    filter(!is.na(mantel_r)) %>%
    transmute(compound = clean_compound_name(name),
              r    = sprintf("%.3f", mantel_r),
              Padj = sprintf("%.3f", p_adj))
  if (!is.null(out_file)) readr::write_tsv(tab, out_file)
  tab
}

# ---- Tip metadata for the bacterial tree ----------------------------------
#
# Phylum and MAG/reference status for every tip, from the
# "<phylum>||<MAG name>" headers of the unfiltered FBA table (so all 498
# bacterial MAGs keep their phylum, including the one without a model).
build_tip_metadata <- function(tree_file, fba_file) {
  tips <- read_tree(tree_file)$tip.label
  cols <- genome_columns(readr::read_tsv(fba_file, n_max = 0, show_col_types = FALSE))
  mag_meta <- tibble::tibble(phylum = sub("\\|\\|.*$", "", cols), label = mag_names(cols))

  tibble::tibble(label = tips) %>%
    left_join(mag_meta, by = "label") %>%
    mutate(source        = if_else(is.na(phylum), "Reference", "MAG"),
           display_label = if_else(source == "Reference", label, NA_character_))
}

# ---- Figure 1: tree with GC, genome-size, and flux rings ------------------
#
# Rings (inner to outer): GC content, log10(genome size), then one ring per
# compound in `compound_ids`, in that order. Rooted on Patescibacteriota;
# internal nodes with ultrafast bootstrap >= 90 marked with a black dot.
plot_tree_with_rings <- function(tree_file, tip_metadata, fba, compound_ids,
                                 ring_width = 0.08, ring_offset = 0.05) {
  tree      <- read_tree(tree_file)
  seed_tips <- tip_metadata$label[tip_metadata$phylum %in% PARAMS$root_phylum_bacteria]
  outgroup  <- ape::extract.clade(tree, ape::getMRCA(tree, seed_tips))$tip.label
  tree      <- ape::root(tree, outgroup = outgroup, resolve.root = TRUE)

  phyla <- sort(unique(na.omit(tip_metadata$phylum)))
  p <- ggtree(tree, layout = "fan", open.angle = 10, size = 0.6, aes(color = phylum)) %<+% tip_metadata +
    scale_color_manual(values = setNames(hue_pal()(length(phyla)), phyla), breaks = phyla,
                       na.value = "grey70", name = "Phylum") +
    geom_tiplab(aes(label = display_label), size = 1.0, offset = 0.02) +
    geom_treescale(color = "black") +
    theme(legend.position = "right")

  support  <- suppressWarnings(as.numeric(p$data$label))
  bs_nodes <- p$data[!p$data$isTip & !is.na(support) & support >= 90, ]

  checkm2 <- readr::read_tsv(INPUT$checkm2, show_col_types = FALSE) %>%
    filter(Name %in% tree$tip.label)
  gc_dat <- transmute(checkm2, label = Name, ring_x = "value", value = GC_Content)
  gs_dat <- transmute(checkm2, label = Name, ring_x = "value", value = log10(Genome_Size))

  quality_layers <- list(
    geom_fruit(data = gc_dat, geom = geom_tile, mapping = aes(y = label, x = ring_x, fill = value),
               pwidth = ring_width, offset = ring_offset, axis.params = list(axis = "none")),
    scale_fill_gradient(low = "grey95", high = "steelblue4", name = "GC content"),
    new_scale_fill(),
    geom_fruit(data = gs_dat, geom = geom_tile, mapping = aes(y = label, x = ring_x, fill = value),
               pwidth = ring_width, offset = ring_offset, axis.params = list(axis = "none")),
    scale_fill_gradient(low = "grey95", high = "forestgreen", name = "log10(Genome size)"),
    new_scale_fill()
  )

  flux <- normalized_flux_matrix(fba)[compound_ids, , drop = FALSE]
  flux_layers <- lapply(rownames(flux), function(id) {
    dat <- tibble::tibble(label = colnames(flux), ring_x = "value", value = flux[id, ])
    geom_fruit(data = dat, geom = geom_tile, mapping = aes(y = label, x = ring_x, fill = value),
               pwidth = ring_width, offset = ring_offset, axis.params = list(axis = "none"))
  })

  ring_names <- c("GC content", "log10(Genome size)",
                  clean_compound_name(fba$name[match(compound_ids, fba$compound_id)]))

  Reduce(`+`, c(quality_layers, flux_layers), init = p) +
    scale_fill_gradient2(low = "purple4", mid = "grey85", high = "darkorange",
                         midpoint = 0, limits = c(-1, 1), name = "Normalized\nflux") +
    geom_point(data = bs_nodes, aes(x = x, y = y), inherit.aes = FALSE,
               shape = 16, size = 1.2, color = "black") +
    labs(caption = paste0("Rings (inner to outer): ",
                          paste(seq_along(ring_names), ring_names, sep = " = ", collapse = " | "))) +
    theme(plot.caption = element_text(hjust = 0, size = 8))
}

# ---- Figure 2 / Supplementary Figure 1: Mantel scatter panels --------------
#
# One panel per compound: patristic distance vs. |difference| in normalized
# flux for every MAG pair, as 2D bin densities (40 x 40, log10 count) with an
# OLS line and 95% band for visual reference. Panels labeled with Mantel r
# and P_adj and ordered by decreasing r.
plot_mantel_scatter <- function(tree_file, fba, mantel, compound_ids, ncol) {
  flux    <- normalized_flux_matrix(fba)
  D_phylo <- mag_patristic_distances(tree_file, colnames(flux))
  upper   <- upper.tri(D_phylo)

  dat <- dplyr::bind_rows(lapply(compound_ids, function(id) {
    D_compound <- as.matrix(stats::dist(flux[id, ]))
    data.frame(compound_id = id, D_phylo = D_phylo[upper], D_compound = D_compound[upper])
  })) %>%
    left_join(mantel, by = "compound_id") %>%
    mutate(panel = sprintf("%s\nr = %.3f (Padj = %.3f)", clean_compound_name(name), mantel_r, p_adj),
           panel = factor(panel, levels = unique(panel[order(-mantel_r)])))

  ggplot(dat, aes(x = D_phylo, y = D_compound)) +
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

# ---- Exchange patterns and GC content -------------------------------------
#
# Logical vectors (named by MAG) for the patterns discussed in the main text.
exchange_patterns <- function(fba) {
  cols <- genome_columns(fba)
  flux <- function(compound) setNames(unlist(fba[fba$name == compound, cols]), mag_names(cols))
  bcaa_import <- flux("L-Leucine [e0]") > 0 | flux("L-Isoleucine [e0]") > 0
  pyr_eth_sec <- flux("Pyruvate [e0]") < 0 & flux("Ethanol [e0]") < 0
  list(
    bcaa_import      = bcaa_import,
    bcaa_secretion   = flux("L-Leucine [e0]") < 0 | flux("L-Isoleucine [e0]") < 0,
    pyr_eth_secretion = pyr_eth_sec,
    acetaldehyde_uptake = flux("Acetaldehyde [e0]") > 0,
    combined         = bcaa_import & pyr_eth_sec,
    phylum           = setNames(sub("^p__", "", sub("\\|\\|.*$", "", cols)), mag_names(cols))
  )
}

# Count and fraction of MAGs with `pattern` per phylum (phyla with >= min_n
# MAGs), highest fraction first.
pattern_by_phylum <- function(pattern, phylum, min_n = 10) {
  tibble::tibble(phylum = phylum, hit = pattern) %>%
    group_by(phylum) %>%
    summarise(n_mags = n(), n_pattern = sum(hit), fraction = n_pattern / n_mags, .groups = "drop") %>%
    filter(n_mags >= min_n) %>%
    arrange(desc(fraction))
}

# Two-sided Wilcoxon rank-sum (Mann-Whitney U) test of CheckM2 GC content
# between MAGs with `pattern` and all other bacterial MAGs.
gc_content_test <- function(pattern, phylum, out_file = NULL) {
  checkm2 <- readr::read_tsv(INPUT$checkm2, show_col_types = FALSE)
  gc <- setNames(checkm2$GC_Content, checkm2$Name)[names(pattern)]
  test <- wilcox.test(gc[pattern], gc[!pattern], alternative = "two.sided", exact = FALSE)

  result <- tibble::tibble(
    n_pattern = sum(pattern), n_other = sum(!pattern),
    gc_min_pattern = min(gc[pattern]), gc_max_pattern = max(gc[pattern]),
    gc_median_pattern = median(gc[pattern]), gc_median_other = median(gc[!pattern]),
    W = unname(test$statistic), p_value = test$p.value
  )
  by_phylum <- tibble::tibble(phylum = phylum[pattern], gc = gc[pattern]) %>%
    group_by(phylum) %>%
    summarise(n = n(), gc_median = median(gc), .groups = "drop") %>%
    arrange(desc(n))

  if (!is.null(out_file)) readr::write_tsv(result, out_file)
  list(test = result, by_phylum = by_phylum)
}
