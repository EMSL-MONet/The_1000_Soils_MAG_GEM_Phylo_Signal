# ============================================================
# 03_build_tip_metadata.R
#
# Per-tip metadata for the bacterial tree (Figure 1): phylum and source
# (MAG vs. GTDB reference) for every tip, parsed from the
# "<phylum>||<MAG name>" headers of the compiled FBA table.
#
# Uses the unfiltered table so that all 498 bacterial MAGs keep their
# phylum, including the one without an FBA model (it is still a tip in
# the tree; it simply has no tiles in the flux rings).
#
# Output columns: label, phylum, source, display_label
# (display_label = tip label for references, NA for MAGs).
# ============================================================

source(file.path("R", "config.R"))
source(file.path("R", "utils.R"))

tips <- read_tree(INPUT$tree_bacteria)$tip.label
cols <- genome_columns(readr::read_tsv(INPUT$fba_bacteria, n_max = 0, show_col_types = FALSE))

mag_meta <- tibble::tibble(phylum = sub("\\|\\|.*$", "", cols), label = mag_names(cols))

missing <- setdiff(mag_meta$label, tips)
if (length(missing) > 0) stop("MAGs not found as tree tips: ", paste(missing, collapse = ", "))

tip_metadata <- tibble::tibble(label = tips) %>%
  left_join(mag_meta, by = "label") %>%
  mutate(source        = if_else(is.na(phylum), "Reference", "MAG"),
         display_label = if_else(source == "Reference", label, NA_character_))

message(sprintf("%d tips (%d MAG, %d Reference)", nrow(tip_metadata),
                sum(tip_metadata$source == "MAG"), sum(tip_metadata$source == "Reference")))
readr::write_tsv(tip_metadata, OUTPUT$tip_metadata_bacteria)
