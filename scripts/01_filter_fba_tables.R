# ============================================================
# 01_filter_fba_tables.R
#
# Prepares the compiled KBase FBA exchange-flux tables for analysis.
#
#   1. Drops MAGs with no FBA output at all (every compound NA). For
#      bacteria this removes PPRH-TOP_bin.5.o, for which no genome-scale
#      model could be generated (550 of 551 MAGs have models).
#   2. Drops compounds missing (NA) in >= one-third of the remaining MAGs
#      (threshold = ceiling(0.333 x number of MAGs)).
#   3. Sets the remaining NAs to zero (compound present in the MAG
#      collection but not exchanged by that MAG's model).
#   4. Sorts compounds by mean flux across MAGs, highest first; ties
#      (e.g., trace ions with identical profiles) keep compound-ID order.
#
# Bacteria: 254 -> 50 compounds across 497 MAGs.
# Archaea:  141 -> 67 compounds across 53 MAGs.
# ============================================================

source(file.path("R", "config.R"))
source(file.path("R", "utils.R"))

filter_fba_table <- function(in_file, out_file, domain) {
  fba  <- readr::read_tsv(in_file, show_col_types = FALSE)
  cols <- genome_columns(fba)

  no_model <- cols[colSums(!is.na(fba[, cols])) == 0]
  if (length(no_model) > 0) {
    message(domain, ": dropping MAG(s) with no FBA output: ",
            paste(mag_names(no_model), collapse = ", "))
    cols <- setdiff(cols, no_model)
  }

  threshold <- ceiling(length(cols) * PARAMS$missing_fraction_cutoff)
  n_missing <- rowSums(is.na(fba[, cols]))
  kept <- fba[n_missing < threshold, c("compound_id", "name", cols)]
  kept[cols] <- lapply(kept[cols], function(x) replace(x, is.na(x), 0))
  kept$avg_uptake <- rowMeans(kept[, cols])
  kept <- kept[order(-kept$avg_uptake), ]

  message(sprintf("%s: %d MAGs; %d of %d compounds retained (NA threshold = %d)",
                  domain, length(cols), nrow(kept), nrow(fba), threshold))
  readr::write_tsv(kept, out_file)
  invisible(kept)
}

filter_fba_table(INPUT$fba_bacteria, OUTPUT$fba_filtered_bacteria, "Bacteria")
filter_fba_table(INPUT$fba_archaea,  OUTPUT$fba_filtered_archaea,  "Archaea")
