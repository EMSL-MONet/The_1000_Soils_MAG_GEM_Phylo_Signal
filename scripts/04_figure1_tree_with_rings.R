# ============================================================
# 04_figure1_tree_with_rings.R
#
# Figure 1: circular bacterial phylogeny (MAGs + GTDB references) with
# concentric rings, innermost first:
#   - GC content and log10(genome size) from CheckM2
#   - normalized flux of each compound with a significant whole-tree
#     Mantel signal (P_adj < 0.05), ordered by Mantel r, highest first
#
# Flux tiles use one diverging scale shared by all rings (purple =
# secretion, grey = no net exchange, orange = uptake; limits -1 to 1),
# after per-compound symmetric min-max scaling. Reference genomes and
# the MAG without an FBA model have no flux tiles.
#
# Tree: rooted on Patescibacteriota; internal nodes with ultrafast
# bootstrap >= 90 marked with a black dot; reference tips labeled.
#
# Output is the raw ggtree figure; final layout (legend arrangement,
# ring numbering) was finished in Adobe Illustrator.
# ============================================================

source(file.path("R", "config.R"))
source(file.path("R", "utils.R"))
suppressPackageStartupMessages({
  library(ggplot2)
  library(ggtree)
  library(ggtreeExtra)
  library(ggnewscale)
  library(tidyr)
  library(scales)
})

ring_width  <- 0.08
ring_offset <- 0.05

# ---- Compounds to plot ----
mantel <- readr::read_tsv(OUTPUT$mantel_bacteria, show_col_types = FALSE)
selected <- mantel %>%
  filter(!is.na(p_adj), p_adj < PARAMS$padj_cutoff) %>%
  arrange(desc(mantel_r))
message("Rings: ", paste(clean_compound_name(selected$name), collapse = ", "))

# ---- Tree, rooted on the outgroup phylum ----
tree     <- read_tree(INPUT$tree_bacteria)
metadata <- readr::read_tsv(OUTPUT$tip_metadata_bacteria, show_col_types = FALSE)

seed_tips <- metadata$label[metadata$phylum %in% PARAMS$root_phylum_bacteria]
outgroup  <- ape::extract.clade(tree, ape::getMRCA(tree, seed_tips))$tip.label
tree      <- ape::root(tree, outgroup = outgroup, resolve.root = TRUE)

phyla <- sort(unique(na.omit(metadata$phylum)))
p <- ggtree(tree, layout = "fan", open.angle = 10, size = 0.6, aes(color = phylum)) %<+% metadata +
  scale_color_manual(values = setNames(hue_pal()(length(phyla)), phyla), breaks = phyla,
                     na.value = "grey70", name = "Phylum") +
  geom_tiplab(aes(label = display_label), size = 1.0, offset = 0.02) +
  geom_treescale(color = "black") +
  theme(legend.position = "right")

# Ultrafast bootstrap support is stored as the internal-node label.
support  <- suppressWarnings(as.numeric(p$data$label))
bs_nodes <- p$data[!p$data$isTip & !is.na(support) & support >= 90, ]

# ---- Genome-quality rings ----
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

# ---- Flux rings ----
fba  <- readr::read_tsv(OUTPUT$fba_filtered_bacteria, show_col_types = FALSE)
flux <- normalized_flux_matrix(fba)[selected$compound_id, , drop = FALSE]

flux_layers <- lapply(rownames(flux), function(id) {
  dat <- tibble::tibble(label = colnames(flux), ring_x = "value", value = flux[id, ])
  geom_fruit(data = dat, geom = geom_tile, mapping = aes(y = label, x = ring_x, fill = value),
             pwidth = ring_width, offset = ring_offset, axis.params = list(axis = "none"))
})

ring_names <- c("GC content", "log10(Genome size)", clean_compound_name(selected$name))

fig <- Reduce(`+`, c(quality_layers, flux_layers), init = p) +
  scale_fill_gradient2(low = "purple4", mid = "grey85", high = "darkorange",
                       midpoint = 0, limits = c(-1, 1), name = "Normalized\nuptake") +
  geom_point(data = bs_nodes, aes(x = x, y = y), inherit.aes = FALSE,
             shape = 16, size = 1.2, color = "black") +
  labs(caption = paste0("Rings (inner to outer): ",
                        paste(seq_along(ring_names), ring_names, sep = " = ", collapse = " | "))) +
  theme(plot.caption = element_text(hjust = 0, size = 8))

out <- file.path(FIG_DIR, "Figure_1_bacteria_tree_with_rings")
ggsave(paste0(out, ".pdf"), fig, width = 28, height = 28, limitsize = FALSE)
ggsave(paste0(out, ".png"), fig, width = 28, height = 28, dpi = 300, limitsize = FALSE)
