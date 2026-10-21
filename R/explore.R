# ============================================================
# explore.R
#
# Helpers for exploring the filtered FBA exchange fluxes
# (Part 2 of analysis.Rmd). Fluxes are in mmol gDW^-1 h^-1;
# positive = uptake, negative = secretion, 0 = no net exchange.
# ============================================================

# Long table: one row per MAG x compound, with phylum.
fba_long <- function(fba) {
  cols <- genome_columns(fba)
  fba %>%
    select(compound_id, name, all_of(cols)) %>%
    pivot_longer(all_of(cols), names_to = "column", values_to = "flux") %>%
    mutate(phylum   = sub("^p__", "", sub("\\|\\|.*$", "", column)),
           mag      = mag_names(column),
           compound = clean_compound_name(name),
           exchange = case_when(flux > 0 ~ "uptake", flux < 0 ~ "secretion", TRUE ~ "none")) %>%
    select(mag, phylum, compound_id, compound, flux, exchange)
}

# The compounds a single MAG takes up and secretes, largest fluxes first.
mag_exchange_profile <- function(long, mag_name) {
  hits <- long %>% filter(mag == mag_name)
  if (nrow(hits) == 0) stop("MAG not found: ", mag_name)
  hits %>%
    filter(exchange != "none") %>%
    arrange(exchange, desc(abs(flux))) %>%
    select(exchange, compound, flux)
}

# For one compound: per phylum, how many MAGs take it up, secrete it, or
# neither, and the median flux. Phyla with fewer than `min_n` MAGs are
# pooled as "Other".
metabolite_by_phylum <- function(long, compound_name, min_n = 10) {
  hits <- long %>% filter(compound == compound_name)
  if (nrow(hits) == 0) stop("Compound not found: ", compound_name)
  sizes <- hits %>% count(phylum)
  hits %>%
    mutate(phylum = if_else(phylum %in% sizes$phylum[sizes$n >= min_n], phylum, "Other")) %>%
    group_by(phylum) %>%
    summarise(n_mags      = n(),
              uptake      = sum(exchange == "uptake"),
              secretion   = sum(exchange == "secretion"),
              none        = sum(exchange == "none"),
              median_flux = median(flux),
              .groups = "drop") %>%
    arrange(phylum == "Other", desc(n_mags))
}

plot_metabolite_by_phylum <- function(summary_tbl, compound_name) {
  summary_tbl %>%
    pivot_longer(c(uptake, secretion, none), names_to = "exchange", values_to = "n") %>%
    mutate(exchange = factor(exchange, levels = c("uptake", "none", "secretion")),
           phylum   = factor(phylum, levels = rev(summary_tbl$phylum))) %>%
    ggplot(aes(x = n, y = phylum, fill = exchange)) +
    geom_col(position = "fill") +
    scale_fill_manual(values = c(uptake = "darkorange", none = "grey80", secretion = "purple4")) +
    scale_x_continuous(labels = scales::percent) +
    labs(x = "Fraction of MAGs", y = NULL, fill = NULL,
         title = paste("Predicted exchange of", compound_name, "by phylum")) +
    theme_minimal(base_size = 11)
}
